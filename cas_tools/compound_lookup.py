"""Web-backed compound metadata lookup for validated CAS registry numbers."""

from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, Mapping, Optional, Tuple
from urllib.parse import quote
from urllib.request import Request, urlopen

from .cas_number_extractor import is_valid_cas_number

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
NCI_RESOLVER_BASE_URL = "https://cactus.nci.nih.gov/chemical/structure"
DEFAULT_TIMEOUT_SECONDS = 10.0
DEFAULT_MAX_SYNONYMS = 12

FetchBytes = Callable[[str, float], bytes]


@dataclass(frozen=True)
class CompoundLookupResult:
    """Partial or complete metadata resolved from a CAS number."""

    cas: str
    status: str
    canonical_name: Optional[str] = None
    iupac_name: Optional[str] = None
    abbreviation: Optional[str] = None
    smiles: Optional[str] = None
    formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    inchi_key: Optional[str] = None
    pubchem_cid: Optional[int] = None
    density: Optional[float] = None
    boiling_point_c: Optional[float] = None
    melting_point_c: Optional[float] = None
    substance_kind: Optional[str] = None
    synonyms: Tuple[str, ...] = ()
    source_ids: Tuple[str, ...] = ()
    source_urls: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()

    @property
    def found(self) -> bool:
        return self.status in {"resolved", "partial"}


def _default_fetch_bytes(url: str, timeout: float) -> bytes:
    request = Request(
        url,
        headers={
            "Accept": "application/json, text/plain;q=0.9",
            "User-Agent": "ConditionAgent-CAS-Lookup/1.0",
        },
    )
    with urlopen(request, timeout=timeout) as response:
        return response.read()


def _safe_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _walk_sections(value: Any) -> Iterable[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        if "TOCHeading" in value:
            yield value
        for nested in value.values():
            yield from _walk_sections(nested)
    elif isinstance(value, list):
        for nested in value:
            yield from _walk_sections(nested)


def _section_information(payload: Mapping[str, Any], heading: str) -> Tuple[Mapping[str, Any], ...]:
    for section in _walk_sections(payload):
        if str(section.get("TOCHeading") or "") == heading:
            return tuple(
                item
                for item in section.get("Information") or ()
                if isinstance(item, Mapping)
            )
    return ()


def _information_strings(information: Mapping[str, Any]) -> Tuple[str, ...]:
    value = information.get("Value") or {}
    if not isinstance(value, Mapping):
        return ()
    strings = []
    direct = value.get("String")
    if direct:
        strings.append(str(direct))
    for item in value.get("StringWithMarkup") or ():
        if isinstance(item, Mapping) and item.get("String"):
            strings.append(str(item["String"]))
    return tuple(strings)


def _temperature_c(payload: Mapping[str, Any], heading: str) -> Optional[float]:
    information = _section_information(payload, heading)
    for item in information:
        value = item.get("Value") or {}
        if not isinstance(value, Mapping):
            continue
        unit = str(value.get("Unit") or "")
        numbers = value.get("Number") or ()
        if numbers and ("°C" in unit or unit.casefold() in {"c", "deg c"}):
            parsed = _safe_float(numbers[0])
            if parsed is not None:
                return round(parsed, 6)
    celsius_pattern = re.compile(r"(-?\d+(?:\.\d+)?)\s*(?:°|deg\s*)?C\b", re.I)
    fahrenheit_pattern = re.compile(r"(-?\d+(?:\.\d+)?)\s*(?:°|deg\s*)?F\b", re.I)
    for item in information:
        for text in _information_strings(item):
            match = celsius_pattern.search(text)
            if match:
                return round(float(match.group(1)), 6)
    for item in information:
        for text in _information_strings(item):
            match = fahrenheit_pattern.search(text)
            if match:
                return round((float(match.group(1)) - 32.0) * 5.0 / 9.0, 6)
    return None


def _density(payload: Mapping[str, Any]) -> Optional[float]:
    information = _section_information(payload, "Density")
    unit_pattern = re.compile(
        r"(-?\d+(?:\.\d+)?)\s*(?:g\s*/\s*(?:mL|cm3|cu\s*cm)|g\.cm-3)",
        re.I,
    )
    relative_pattern = re.compile(
        r"(?:relative\s+density|specific\s+gravity)[^\d<]*<?\s*(\d+(?:\.\d+)?)",
        re.I,
    )
    for item in information:
        for text in _information_strings(item):
            match = unit_pattern.search(text)
            if match:
                return round(float(match.group(1)), 6)
    for item in information:
        for text in _information_strings(item):
            match = relative_pattern.search(text)
            if match:
                return round(float(match.group(1)), 6)
    for item in information:
        for text in _information_strings(item):
            match = re.match(r"\s*(\d+(?:\.\d+)?)\b", text)
            if match:
                candidate = float(match.group(1))
                if 0.0 < candidate < 30.0:
                    return round(candidate, 6)
    return None


def _substance_kind(payload: Mapping[str, Any]) -> Optional[str]:
    for item in _section_information(payload, "Physical Description"):
        for text in _information_strings(item):
            lowered = text.casefold()
            if "liquid" in lowered:
                return "liquid"
            if "powder" in lowered:
                return "powder"
            if "solid" in lowered:
                return "solid"
            if "gas" in lowered or "vapor" in lowered:
                return "gas"
    return None


def _looks_like_registry_identifier(value: str) -> bool:
    text = value.strip()
    if is_valid_cas_number(text):
        return True
    if re.fullmatch(r"(?:CID|SID|CHEBI|DTXSID|UNII|NSC)[:\s-]*[A-Z0-9-]+", text, re.I):
        return True
    return bool(re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", text))


def _select_synonyms(
    values: Iterable[Any],
    *,
    excluded: Iterable[str],
    maximum: int,
) -> Tuple[str, ...]:
    excluded_keys = {str(value).strip().casefold() for value in excluded if value}
    selected = []
    seen = set(excluded_keys)
    for raw in values:
        value = " ".join(str(raw or "").split())
        key = value.casefold()
        if not value or key in seen or _looks_like_registry_identifier(value):
            continue
        if len(value) > 120:
            continue
        seen.add(key)
        selected.append(value)
        if len(selected) >= maximum:
            break
    return tuple(selected)


def _abbreviation(values: Iterable[str]) -> Optional[str]:
    for value in values:
        text = value.strip()
        if not 2 <= len(text) <= 12 or any(character.isspace() for character in text):
            continue
        if not re.search(r"[A-Za-z]", text):
            continue
        if re.search(r"[A-Z0-9]", text[1:]) or (text.isupper() and len(text) <= 6):
            return text
    return None


class CompoundLookupClient:
    """Small injectable HTTP client for PubChem with an NCI fallback."""

    def __init__(
        self,
        *,
        fetch_bytes: FetchBytes = _default_fetch_bytes,
        timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
        max_synonyms: int = DEFAULT_MAX_SYNONYMS,
    ) -> None:
        self._fetch_bytes = fetch_bytes
        self.timeout_seconds = float(timeout_seconds)
        self.max_synonyms = int(max_synonyms)

    def _json(self, url: str) -> Mapping[str, Any]:
        content = self._fetch_bytes(url, self.timeout_seconds)
        payload = json.loads(content.decode("utf-8"))
        if not isinstance(payload, Mapping):
            raise ValueError("Lookup response must be a JSON object")
        return payload

    def _text(self, url: str) -> str:
        return self._fetch_bytes(url, self.timeout_seconds).decode("utf-8").strip()

    def _parallel_json(
        self,
        urls: Mapping[str, str],
    ) -> Tuple[Dict[str, Mapping[str, Any]], Dict[str, Exception]]:
        values: Dict[str, Mapping[str, Any]] = {}
        errors: Dict[str, Exception] = {}
        with ThreadPoolExecutor(max_workers=max(1, len(urls))) as executor:
            futures = {
                executor.submit(self._json, url): key
                for key, url in urls.items()
            }
            for future in as_completed(futures):
                key = futures[future]
                try:
                    values[key] = future.result()
                except Exception as error:
                    errors[key] = error
        return values, errors

    def _parallel_text(
        self,
        urls: Mapping[str, str],
    ) -> Tuple[Dict[str, str], Dict[str, Exception]]:
        values: Dict[str, str] = {}
        errors: Dict[str, Exception] = {}
        with ThreadPoolExecutor(max_workers=max(1, len(urls))) as executor:
            futures = {
                executor.submit(self._text, url): key
                for key, url in urls.items()
            }
            for future in as_completed(futures):
                key = futures[future]
                try:
                    values[key] = future.result()
                except Exception as error:
                    errors[key] = error
        return values, errors

    def _pubchem_view(self, cid: int, heading: str) -> Mapping[str, Any]:
        url = (
            f"{PUBCHEM_VIEW_BASE_URL}/data/compound/{cid}/JSON"
            f"?heading={quote(heading)}"
        )
        return self._json(url)

    def lookup(self, cas: str) -> CompoundLookupResult:
        """Look up as many supported compound fields as available."""
        cas = str(cas or "").strip()
        if not is_valid_cas_number(cas):
            return CompoundLookupResult(
                cas=cas,
                status="invalid_identifier",
                warnings=("INVALID_CAS",),
            )
        encoded_cas = quote(cas, safe="")
        property_url = (
            f"{PUBCHEM_BASE_URL}/compound/name/{encoded_cas}/property/"
            "Title,IUPACName,MolecularFormula,MolecularWeight,SMILES,"
            "ConnectivitySMILES,InChIKey/JSON"
        )
        synonym_url = (
            f"{PUBCHEM_BASE_URL}/compound/name/{encoded_cas}/synonyms/JSON"
        )
        warnings = []
        source_urls = []
        source_ids = []
        properties: Dict[str, Any] = {}
        raw_synonyms: Tuple[Any, ...] = ()
        core_payloads, core_errors = self._parallel_json(
            {"properties": property_url, "synonyms": synonym_url}
        )
        payload = core_payloads.get("properties")
        if payload is not None:
            property_values = (
                payload.get("PropertyTable", {}).get("Properties", ())
                if isinstance(payload.get("PropertyTable"), Mapping)
                else ()
            )
            if property_values:
                properties = dict(property_values[0])
                if len(property_values) > 1:
                    warnings.append("MULTIPLE_PUBCHEM_COMPOUNDS_FIRST_SELECTED")
                source_urls.append(property_url)
        elif "properties" in core_errors:
            error = core_errors["properties"]
            warnings.append(f"PUBCHEM_PROPERTY_LOOKUP_FAILED:{type(error).__name__}")
        synonym_payload = core_payloads.get("synonyms")
        if synonym_payload is not None:
            information = synonym_payload.get("InformationList", {}).get("Information", ())
            if information:
                raw_synonyms = tuple(information[0].get("Synonym") or ())
                source_urls.append(synonym_url)
        elif "synonyms" in core_errors:
            error = core_errors["synonyms"]
            warnings.append(f"PUBCHEM_SYNONYM_LOOKUP_FAILED:{type(error).__name__}")

        cid_value = properties.get("CID")
        cid = int(cid_value) if isinstance(cid_value, (int, float)) else None
        canonical_name = str(properties.get("Title") or "").strip() or None
        iupac_name = str(properties.get("IUPACName") or "").strip() or None
        smiles = str(
            properties.get("SMILES") or properties.get("ConnectivitySMILES") or ""
        ).strip() or None
        formula = str(properties.get("MolecularFormula") or "").strip() or None
        molecular_weight = _safe_float(properties.get("MolecularWeight"))
        inchi_key = str(properties.get("InChIKey") or "").strip() or None
        if cid is not None:
            source_ids.append(f"pubchem:cid:{cid}")
            source_urls.append(f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}")

        density = None
        boiling_point_c = None
        melting_point_c = None
        substance_kind = None
        if cid is not None:
            headings = (
                "Density",
                "Boiling Point",
                "Melting Point",
                "Physical Description",
            )
            view_urls = {
                heading: (
                    f"{PUBCHEM_VIEW_BASE_URL}/data/compound/{cid}/JSON"
                    f"?heading={quote(heading)}"
                )
                for heading in headings
            }
            view_payloads, view_errors = self._parallel_json(view_urls)
            for heading in headings:
                view_payload = view_payloads.get(heading)
                if view_payload is None:
                    error = view_errors.get(heading)
                    warnings.append(
                        f"PUBCHEM_{heading.upper().replace(' ', '_')}_LOOKUP_FAILED:"
                        f"{type(error).__name__ if error is not None else 'UnknownError'}"
                    )
                    continue
                if heading == "Density":
                    density = _density(view_payload)
                elif heading == "Boiling Point":
                    boiling_point_c = _temperature_c(view_payload, heading)
                elif heading == "Melting Point":
                    melting_point_c = _temperature_c(view_payload, heading)
                else:
                    substance_kind = _substance_kind(view_payload)
                source_urls.append(view_urls[heading])
            if any(
                value is not None
                for value in (
                    density,
                    boiling_point_c,
                    melting_point_c,
                    substance_kind,
                )
            ):
                warnings.append("PUBCHEM_EXPERIMENTAL_FIELDS_REQUIRE_REVIEW")

        if not any((canonical_name, smiles, formula, inchi_key)):
            fallback = self._nci_fallback(cas)
            canonical_name = (
                canonical_name
                or fallback.get("canonical_name")
                or fallback.get("iupac_name")
            )
            iupac_name = iupac_name or fallback.get("iupac_name")
            smiles = smiles or fallback.get("smiles")
            formula = formula or fallback.get("formula")
            molecular_weight = molecular_weight or _safe_float(
                fallback.get("molecular_weight")
            )
            inchi_key = inchi_key or fallback.get("inchi_key")
            fallback_synonyms = fallback.get("synonyms") or ()
            raw_synonyms = (*raw_synonyms, *fallback_synonyms)
            source_urls.extend(fallback.get("source_urls") or ())
            if fallback:
                source_ids.append("nci:cactus")
            else:
                warnings.append("NCI_FALLBACK_NOT_FOUND")

        synonyms = _select_synonyms(
            raw_synonyms,
            excluded=(cas, canonical_name or "", iupac_name or ""),
            maximum=self.max_synonyms,
        )
        abbreviation = _abbreviation(str(value) for value in raw_synonyms)
        status = "resolved" if all(
            (canonical_name, smiles, formula, molecular_weight)
        ) else "partial"
        if not any((canonical_name, smiles, formula, inchi_key, synonyms)):
            status = "not_found"
        return CompoundLookupResult(
            cas=cas,
            status=status,
            canonical_name=canonical_name,
            iupac_name=iupac_name,
            abbreviation=abbreviation,
            smiles=smiles,
            formula=formula,
            molecular_weight=molecular_weight,
            inchi_key=inchi_key,
            pubchem_cid=cid,
            density=density,
            boiling_point_c=boiling_point_c,
            melting_point_c=melting_point_c,
            substance_kind=substance_kind,
            synonyms=synonyms,
            source_ids=tuple(dict.fromkeys(source_ids)),
            source_urls=tuple(dict.fromkeys(source_urls)),
            warnings=tuple(dict.fromkeys(warnings)),
        )

    def _nci_fallback(self, cas: str) -> Dict[str, Any]:
        encoded = quote(cas, safe="")
        outputs = {
            "iupac_name": "iupac_name",
            "smiles": "smiles",
            "formula": "formula",
            "molecular_weight": "mw",
            "inchi_key": "stdinchiKey",
        }
        urls = {
            field: f"{NCI_RESOLVER_BASE_URL}/{encoded}/{representation}"
            for field, representation in outputs.items()
        }
        names_url = f"{NCI_RESOLVER_BASE_URL}/{encoded}/names"
        values, _errors = self._parallel_text({**urls, "names": names_url})
        result: Dict[str, Any] = {"source_urls": []}
        for field, url in urls.items():
            value = values.get(field, "")
            if value and not value.startswith("<"):
                result[field] = value
                result["source_urls"].append(url)
        names = tuple(
            line.strip()
            for line in values.get("names", "").splitlines()
            if line.strip()
        )
        selected_names = _select_synonyms(
            names,
            excluded=(cas,),
            maximum=max(self.max_synonyms + 1, 2),
        )
        if selected_names:
            result["canonical_name"] = selected_names[0]
            result["synonyms"] = selected_names
            result["source_urls"].append(names_url)
        return result if len(result) > 1 else {}


def lookup_compound_by_cas(
    cas: str,
    *,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    max_synonyms: int = DEFAULT_MAX_SYNONYMS,
) -> CompoundLookupResult:
    """Convenience API for one web-backed CAS lookup."""
    return CompoundLookupClient(
        timeout_seconds=timeout_seconds,
        max_synonyms=max_synonyms,
    ).lookup(cas)


__all__ = [
    "DEFAULT_MAX_SYNONYMS",
    "DEFAULT_TIMEOUT_SECONDS",
    "NCI_RESOLVER_BASE_URL",
    "PUBCHEM_BASE_URL",
    "PUBCHEM_VIEW_BASE_URL",
    "CompoundLookupClient",
    "CompoundLookupResult",
    "lookup_compound_by_cas",
]

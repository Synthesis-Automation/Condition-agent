$path = 'chemtools/condition_core.py'
$lines = Get-Content $path -Encoding UTF8
$start = ($lines | Select-String -Pattern '^def _read_dataset_aliases\(' -SimpleMatch | Select-Object -First 1).LineNumber
$end = ($lines | Select-String -Pattern '^\s*\n\s*\n_LIG_BY_CAS' | Select-Object -First 1).LineNumber
$start
$end
$lines[($start-1)..($end-2)] | ForEach-Object { $_ }

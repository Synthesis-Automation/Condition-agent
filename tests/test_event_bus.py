"""Phase 2: EventBus unit tests."""

import pytest
from chem_coworker.event_bus import EventBus, ChemEvent


class TestEventBusSubscribeEmit:
    def test_subscribe_and_emit_roundtrip(self):
        bus = EventBus()
        received = []
        bus.subscribe(ChemEvent.TOOL_DONE, lambda tool_name, elapsed_s, **_: received.append((tool_name, elapsed_s)))
        bus.emit(ChemEvent.TOOL_DONE, tool_name="detect_reaction_type", elapsed_s=0.12)
        assert received == [("detect_reaction_type", 0.12)]

    def test_multiple_handlers_all_fire(self):
        bus = EventBus()
        log = []
        bus.subscribe(ChemEvent.PHASE_START, lambda phase, **_: log.append(f"A:{phase}"))
        bus.subscribe(ChemEvent.PHASE_START, lambda phase, **_: log.append(f"B:{phase}"))
        bus.emit(ChemEvent.PHASE_START, phase="reason")
        assert log == ["A:reason", "B:reason"]

    def test_no_handlers_no_error(self):
        bus = EventBus()
        bus.emit(ChemEvent.TOOL_START, tool_name="foo")  # should not raise

    def test_handler_exception_does_not_abort(self):
        bus = EventBus()
        log = []
        bus.subscribe(ChemEvent.TOOL_DONE, lambda **_: (_ for _ in ()).throw(RuntimeError("boom")))
        bus.subscribe(ChemEvent.TOOL_DONE, lambda tool_name, **_: log.append(tool_name))
        bus.emit(ChemEvent.TOOL_DONE, tool_name="ok_tool", elapsed_s=0.1)
        # Second handler should still have fired
        assert log == ["ok_tool"]

    def test_handler_exception_does_not_abort_v2(self):
        """Simpler exception check without generator tricks."""
        bus = EventBus()
        log = []

        def bad_handler(**kwargs):
            raise ValueError("intentional failure")

        bus.subscribe(ChemEvent.STREAM_TOKEN, bad_handler)
        bus.subscribe(ChemEvent.STREAM_TOKEN, lambda token, **_: log.append(token))
        bus.emit(ChemEvent.STREAM_TOKEN, token="hello")
        assert log == ["hello"]

    def test_handlers_fire_in_subscription_order(self):
        bus = EventBus()
        order = []
        bus.subscribe(ChemEvent.PHASE_END, lambda **_: order.append(1))
        bus.subscribe(ChemEvent.PHASE_END, lambda **_: order.append(2))
        bus.subscribe(ChemEvent.PHASE_END, lambda **_: order.append(3))
        bus.emit(ChemEvent.PHASE_END, phase="synth")
        assert order == [1, 2, 3]


class TestEventBusBoolSemantics:
    def test_empty_bus_is_falsy(self):
        bus = EventBus()
        assert not bus

    def test_bus_with_handler_is_truthy(self):
        bus = EventBus()
        bus.subscribe(ChemEvent.TOOL_DONE, lambda **_: None)
        assert bus

    def test_bus_stays_truthy_after_emit(self):
        bus = EventBus()
        bus.subscribe(ChemEvent.COMPACT_START, lambda **_: None)
        bus.emit(ChemEvent.COMPACT_START)
        assert bus


class TestEventBusDecoratorSyntax:
    def test_on_decorator_registers_handler(self):
        bus = EventBus()
        received = []

        @bus.on(ChemEvent.PRE_SYNTH)
        def handler(hypothesis, confidence, **_):
            received.append((hypothesis, confidence))

        bus.emit(ChemEvent.PRE_SYNTH, hypothesis="Suzuki", confidence=0.9,
                 rationale="test", tools_called=[], plan_revised=False)
        assert received == [("Suzuki", 0.9)]

    def test_on_decorator_returns_original_fn(self):
        bus = EventBus()

        @bus.on(ChemEvent.TOOL_ERROR)
        def my_fn(**_):
            pass

        assert my_fn.__name__ == "my_fn"


class TestEventBusChaining:
    def test_subscribe_returns_self(self):
        bus = EventBus()
        result = bus.subscribe(ChemEvent.TOOL_START, lambda **_: None)
        assert result is bus

    def test_chained_subscribe(self):
        bus = EventBus()
        log = []
        (bus
            .subscribe(ChemEvent.TOOL_DONE, lambda tool_name, **_: log.append(f"done:{tool_name}"))
            .subscribe(ChemEvent.TOOL_ERROR, lambda tool_name, **_: log.append(f"err:{tool_name}"))
        )
        bus.emit(ChemEvent.TOOL_DONE, tool_name="resolve_to_smiles", elapsed_s=0.5)
        bus.emit(ChemEvent.TOOL_ERROR, tool_name="bad_tool", elapsed_s=0.1, error="oops")
        assert log == ["done:resolve_to_smiles", "err:bad_tool"]


class TestChemEventEnum:
    def test_all_expected_events_exist(self):
        expected = {
            "TOOL_START", "TOOL_DONE", "TOOL_ERROR",
            "PHASE_START", "PHASE_END",
            "PRE_SYNTH", "STREAM_TOKEN",
            "COMPACT_START", "COMPACT_END",
        }
        actual = {e.name for e in ChemEvent}
        assert expected == actual

    def test_events_are_distinct(self):
        events = list(ChemEvent)
        assert len(events) == len(set(events))

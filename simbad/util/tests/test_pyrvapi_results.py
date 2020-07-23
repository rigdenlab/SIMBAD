import unittest
from simbad.util import RvapiMetadata


class TestRvapiMetadata(unittest.TestCase):
    def test_add_1(self):
        meta = RvapiMetadata()
        meta.add({})
        for ans, meta in zip([{}], meta.results):
            self.assertDictEqual(ans, meta)

    def test_add_2(self):
        meta = RvapiMetadata()
        meta.add({"key": "board"})
        meta.add({"foo": "bar"})
        for ans, meta in zip([{"key": "board"}, {"foo": "bar"}], meta.results):
            self.assertDictEqual(ans, meta)

    def test_add_3(self):
        meta = RvapiMetadata()
        meta.add({})
        meta.add({"foo": "bar"})
        for ans, meta in zip([{}, {"foo": "bar"}], meta.results):
            self.assertDictEqual(ans, meta)

    def test_add_4(self):
        meta = RvapiMetadata()
        self.assertEqual(0, len(meta.results))

    def test_to_json_1(self):
        meta = RvapiMetadata()
        meta.add({"key": "board"})
        answer = '{"nResults": 1, "results": [{"key": "board"}]}'
        self.assertEqual(answer, meta.to_json())

    def test_to_json_2(self):
        meta = RvapiMetadata()
        answer = '{"nResults": 0, "results": []}'
        self.assertEqual(answer, meta.to_json())

    def test_to_json_3(self):
        meta = RvapiMetadata()
        meta.add({})
        meta.add({"foo": "bar"})
        answer = '{"nResults": 2, "results": [{}, {"foo": "bar"}]}'
        self.assertEqual(answer, meta.to_json())


if __name__ == "__main__":
    unittest.main(verbosity=2)

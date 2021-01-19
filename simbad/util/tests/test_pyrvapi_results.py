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
        answer_py2 = '{"nResults": 1, "results": [{"key": "board"}]}'
        answer_py3 = '{"results": [{"key": "board"}], "nResults": 1}'
        self.assertTrue(meta.to_json() in [answer_py2, answer_py3])

    def test_to_json_2(self):
        meta = RvapiMetadata()
        answer_py2 = '{"nResults": 0, "results": []}'
        answer_py3 = '{"results": [], "nResults": 0}'
        self.assertTrue(meta.to_json() in [answer_py2, answer_py3])

    def test_to_json_3(self):
        meta = RvapiMetadata()
        meta.add({})
        meta.add({"foo": "bar"})
        answer_py2 = '{"nResults": 2, "results": [{}, {"foo": "bar"}]}'
        answer_py3 = '{"results": [{}, {"foo": "bar"}], "nResults": 2}'
        self.assertTrue(meta.to_json() in [answer_py2, answer_py3])


if __name__ == "__main__":
    unittest.main(verbosity=2)

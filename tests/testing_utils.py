# Run all tests on one chosen model, both for efficiency's sake
# and because we'd like to have a single focal model
TEST_MODEL = "model/Rpom_05.xml"

def run_all_tests_in_object(tests_object, test_prefix="test", test_params=None):
    for test_func in [f for f in dir(tests_object) if f.startswith(test_prefix)]:
        test = getattr(tests_object, test_func)
        if test_params is None or test_func not in test_params:
            test()
        else:
            for param_set in test_params[test_func]:
                test(*param_set)

def run_all_tests_in_object(tests_object, test_prefix="test"):
    for test_func in [f for f in dir(tests_object) if f.startswith(test_prefix)]:
        getattr(tests_object, test_func)()

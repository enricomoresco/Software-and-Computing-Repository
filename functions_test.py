import testing_functions as test

test.mass_conservation(100)

test.momentum_conservation(100)

test.verify_BC_evolution(100)

print("#################")
print("ALL TESTS PASSED")
print("#################")
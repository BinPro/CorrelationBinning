#!/usr/bin/env python
import unittest
from nose.tools import assert_almost_equal, assert_equal

from corrbin.classification import false_positive_rate, \
    true_positive_rate

# testing function: false_positive_rate
# False positive = 0, Negative: FP + TN = 2
def test_false_positive_rate1():
    est_positive = [True, True]
    est_negative = [False,False]
    
    result = false_positive_rate(est_positive,est_negative)
    assert_almost_equal(result,0.0)

# False positive = 1, Negative: FP + TN = 3
def test_false_positive_rate2():
    est_positive = [False, True]
    est_negative = [False,False]
    
    result = false_positive_rate(est_positive,est_negative)
    assert_almost_equal(result,1/3.0)

# False positive = 2, Negative: FP + TN = 4
def test_false_positive_rate3():
    est_positive = [False, False]
    est_negative = [False,False]
    
    result = false_positive_rate(est_positive,est_negative)
    assert_almost_equal(result,2/4.0)

# False positive = 2, Negative: FP + TN = 3
def test_false_positive_rate4():
    est_positive = [False, False]
    est_negative = [False,True]
    
    result = false_positive_rate(est_positive,est_negative)
    assert_almost_equal(result,2/3.0)

# False positive = 2, Negative: FP + TN = 2
def test_false_positive_rate5():
    est_positive = [False, False]
    est_negative = [True,True]
    
    result = false_positive_rate(est_positive,est_negative)
    assert_almost_equal(result,2/2.0)

# True positive = 2, Positive: TP + FN = 2
def test_true_positive_rate1():
    est_positive = [True, True]
    est_negative = [False,False]
    
    result = true_positive_rate(est_positive,est_negative)
    assert_almost_equal(result,2/2.0)


# True positive = 2, Positive: TP + FN = 3
def test_true_positive_rate2():
    est_positive = [True, True]
    est_negative = [False,True]
    
    result = true_positive_rate(est_positive,est_negative)
    assert_almost_equal(result,2/3.0)


# True positive = 2, Positive: TP + FN = 4
def test_true_positive_rate3():
    est_positive = [True, True]
    est_negative = [True,True]
    
    result = true_positive_rate(est_positive,est_negative)
    assert_almost_equal(result,2/4.0)


# True positive = 1, Positive: TP + FN = 1
def test_true_positive_rate4():
    est_positive = [False, True]
    est_negative = [False,False]
    
    result = true_positive_rate(est_positive,est_negative)
    assert_almost_equal(result,1/1.0)


# True positive = 0, Positive: TP + FN = 0
def test_true_positive_rate5():
    est_positive = [False,False]
    est_negative = [False,False]
    
    result = true_positive_rate(est_positive,est_negative)
    assert_almost_equal(result,0.0)


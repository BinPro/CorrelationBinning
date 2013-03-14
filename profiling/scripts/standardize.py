#!/usr/bin/env python
import os
from corrbin.score import ExperimentData, RocAxisFun
file_path = os.path.realpath(__file__)
input_file = os.path.abspath(os.path.join(file_path,"..","..","fixtures/score_500000.tsv"))
@profile
def my_fun():
    fun = RocAxisFun("precision").fun
    data = ExperimentData(fun)
    data.load_data_frame(input_file)
    data.standardize()

my_fun()

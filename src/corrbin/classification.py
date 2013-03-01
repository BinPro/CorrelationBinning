import pandas as p
import corrbin.exceptions as exc

# Classifies each contig in dataframe df
# according to cut-off value: q.
# 
# Returns a dataframe with one line for each contig_id
# and a column with contig_id and one with classification.
#
# Clasification is True if max contig score occur
# on the line where contig_genome == compare_genome 
# Output is False if max contig score occur on any
# other line than above.
# If the max contig score is less than q, the classification
# will be negative, otherwise positive.
def classify_bool(d, level):
    df = d.df
    contig_ids = df.contig_id.unique()
    output_real = []
    output_q = []
    for contig_id in contig_ids:
        contig_df = df[df.contig_id == contig_id]
        real_class, max_q_std = classify_contig_bool(contig_df,level)
        output_real.append(real_class)
        output_q.append(max_q_std)
    s_real = p.Series(output_real,index=contig_ids)
    s_max_q = p.Series(output_q, index=contig_ids)
    df = p.DataFrame({"real_classif": s_real,
                      "max_std_p": s_max_q})
    return df


def classify_contig_bool(df,level):
    max_q = df.p_value_standardized.max()
    # Since eval is used later in code, check that level
    # is a safe string.
    if not (level in ["genome","species", "genus","family"]):
        raise exc.LevelError(level)
    real_class = False
    max_rows = df[df.p_value_standardized == max_q]
    if len(max_rows)>1:
        for row in max_rows:
            real_class = \
                (eval("row.contig_" + level) == \
                     eval("row.compare_" + level)) or \
                     real_class
    else:
        real_class = (eval("max_rows.contig_" + level) == \
                          eval("max_rows.compare_" + level))
    return bool(real_class), max_q

# Calculates the data needed to make a ROC-curve
# Input: dataframe with columns "max_std_p" 
# and a boolean as the second column.
# Output: A data frame with two columns, True positive rate
# and false positive rate calculated using each 
# max_std_p value as a threshold respectively.
def calculate_roc(df,axis_funs):
    x_fun = axis_funs.x_fun
    y_fun = axis_funs.y_fun
    # Sort max_p
    sorted_class = [row[1] for row in df.sort(columns="max_std_p", ascending=True).values]
    # Stepwise, calculate the roc_data
    est_positives = []
    y = []
    x = []
    while sorted_class != []:
        est_positives.append(sorted_class.pop())
        y.append(y_fun(est_positives, sorted_class))
        x.append(x_fun(est_positives, sorted_class))
    # collect in a series
    y_s = p.Series(y)
    x_s = p.Series(x)
    df = p.DataFrame({"x":x_s, "y":y_s})
    return df

def true_positive_rate(est_positive, est_negative):
    TP = est_positive.count(True)
    FP = est_positive.count(False)
    TN = est_negative.count(False)
    FN = est_negative.count(True)
    
    if TP == 0:
        return 0.0
    else:
        return float(TP)/(TP + FN)

def false_positive_rate(est_positive, est_negative):
    TP = est_positive.count(True)
    FP = est_positive.count(False)
    TN = est_negative.count(False)
    FN = est_negative.count(True)
    
    if FP == 0:
        return 0.0
    else:
        return float(FP)/(FP+TN)

def precision(est_positive, est_negative):
    TP = est_positive.count(True)
    FP = est_positive.count(False)
    TN = est_negative.count(False)
    FN = est_negative.count(True)
    
    if TP == 0:
        return 0.0
    else:
        return float(TP)/(TP+FP)

def false_discovery_rate(est_positive, est_negative):
    TP = est_positive.count(True)
    FP = est_positive.count(False)
    TN = est_negative.count(False)
    FN = est_negative.count(True)
    
    if FP == 0:
        return 0.0
    else:
        return float(FP)/(TP+FP)

def included_contigs_ratio(est_positive,est_negative):
    total = len(est_positive)+len(est_negative)
    return float(len(est_positive))/total


def accuracy(est_positive,est_negative):
    TP = est_positive.count(True)
    FP = est_positive.count(False)
    TN = est_negative.count(False)
    FN = est_negative.count(True)
    
    total = TP+FP+TN+FN
    return float(TP+TN)/total

def sensitivity(df,q):
    pass

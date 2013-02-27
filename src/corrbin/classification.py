import pandas as p

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
def classify_bool(d):
    df = d.df
    contig_ids = df.contig_id.unique()
    output_real = []
    output_q = []
    for contig_id in contig_ids:
        contig_df = df[df.contig_id == contig_id]
        real_class, max_q_std = classify_contig_bool(contig_df)
        output_real.append(real_class)
        output_q.append(max_q_std)
    s_real = p.Series(output_real,index=contig_ids)
    s_max_q = p.Series(output_q, index=contig_ids)
    df = p.DataFrame({"real_classif": s_real,
                      "max_std_p": s_max_q})
    return df


def classify_contig_bool(df):
    max_q = df.p_value_standardized.max()

    real_class = False
    max_rows = df[df.p_value_standardized == max_q]
    if len(max_rows)>1:
        for row in max_rows:
            real_class = \
                (row.contig_genome == row.compare_genome) or \
                real_class
    else:
        real_class = (max_rows.contig_genome == max_rows.compare_genome)
    return bool(real_class), max_q
        

def sensitivity(df,q):
    pass

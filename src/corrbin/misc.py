class Uniq_id(object):
    def __init__(self, start_val=0, incr=1):
        self.start = start_val
        self.curr_val = start_val
        self.incr = incr

    def id(self):
        val = self.curr_val
        self.curr_val +=self.incr
        return val

def all_but_index(l,i):
    return l[0:i] + l[i+1:]

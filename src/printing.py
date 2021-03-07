def sec_to_str(t):
    # Convert seconds to hours:minutes:seconds
    m, s = divmod(t, 60)
    h, m =  divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)

def _remove_dtype(x):
    # Removes dtype: float64 and dtype: int64 from pandas printouts
    x = str(x)
    x = x.replace('\ndtype: int64', '')
    x = x.replace('\ndtype: float64', '')
    return x

class Logger(object):
    # Lightweight logging.
    def __init__(self, fh):
        self.log_fh = open(fh, 'w')

    def log(self, msg):
        # Print to log file and stdout with a single command.
        self.log_fh.write(msg)
        print(msg)
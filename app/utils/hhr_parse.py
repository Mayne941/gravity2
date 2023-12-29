from io import StringIO
import pandas as pd

def get_template():
    return \
    u"""
    ---|-------------------------------|----|-------|-------|------|-----|----|---------|--------------|
    No Hit                           Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
    """
def hhparse(hhresult_file):
	'''Convert HHpred's text-based output table in to a pandas dataframe'''
	pattern = StringIO(get_template()).readlines()[1]
	colBreaks = [i for i, ch in enumerate(pattern) if ch == '|']
	widths = [j-i for i, j in zip(([0]+colBreaks)[:-1], colBreaks)]
	# df = pd.read_fwf(hhresult_file, skiprows=8,
	#                       nrows=10, header=0, widths=widths)
	df = pd.read_fwf(hhresult_file, colspecs=(3, 30, 5, 7), nrows=10, skiprows=8)
	print(df)
	return df

if __name__ == '__main__':
	hhparse("test.fuck")

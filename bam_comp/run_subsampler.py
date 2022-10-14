from subsampler import CSV, Subsampler


csv_1 = CSV('') # path to first CSV
csv_2 = CSV('') # path to second CSV


subsampler = Subsampler(csv_1, csv_2)
subsampler.run()
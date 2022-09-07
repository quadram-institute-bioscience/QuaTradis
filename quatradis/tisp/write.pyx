def write_plot_file(filename, forward, reverse, p_len):
	with open(filename, 'w') as fh:
		for i in range(0, p_len):
			fh.write(str(int(forward[i]))+' '+str(int(reverse[i]))+"\n")
	return

def write_gzipped_plot_file(filename, forward, reverse, p_len):
	from Bio import bgzf

	with bgzf.BgzfWriter(filename, 'wb') as fh:
		for i in range(0, p_len):
			fh.write(str(int(forward[i]))+' '+str(int(reverse[i]))+"\n")
	return
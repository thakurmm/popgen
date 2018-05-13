# line = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA06984	NA06985	NA06986	NA06989	NA06994	NA07000	NA07037	NA07048	NA07051	NA07056	NA07347	NA07357	NA10847	NA10851	NA11829	NA11830	NA11831	NA11832	NA11840	NA11843	NA11881	NA11892	NA11893	NA11894	NA11918	NA11919	NA11920	NA11930	NA11931	NA11932	NA11933	NA11992	NA11994	NA11995	NA12003	NA12004	NA12005	NA12006	NA12043	NA12044	NA12045	NA12046	NA12058	NA12144	NA12154	NA12155	NA12156	NA12234	NA12249	NA12272	NA12273	NA12275	NA12282	NA12283	NA12286	NA12287	NA12340	NA12341	NA12342	NA12347	NA12348	NA12383	NA12399	NA12400	NA12413	NA12414	NA12489	NA12546	NA12716	NA12717	NA12718	NA12748	NA12749	NA12750	NA12751	NA12760	NA12761	NA12762	NA12763	NA12775	NA12776	NA12777	NA12778	NA12812	NA12813	NA12814	NA12815	NA12827	NA12828	NA12829	NA12830	NA12842	NA12843	NA12872	NA12873	NA12874	NA12878	NA12889	NA12890	NA18486	NA18488	NA18489	NA18498	NA18499	NA18501	NA18502	NA18504	NA18505	NA18507	NA18508	NA18510	NA18511	NA18516	NA18517	NA18519	NA18520	NA18522	NA18523	NA18853	NA18856	NA18858	NA18861	NA18864	NA18865	NA18867	NA18868	NA18870	NA18871	NA18873	NA18874	NA18876	NA18877	NA18878	NA18879	NA18881	NA18907	NA18908	NA18909	NA18910	NA18912	NA18915	NA18916	NA18917	NA18923	NA18924	NA18933	NA18934	NA19092	NA19093	NA19095	NA19096	NA19098	NA19099	NA19102	NA19107	NA19108	NA19113	NA19114	NA19116	NA19117	NA19118	NA19119	NA19121	NA19129	NA19130	NA19131	NA19137	NA19138	NA19141	NA19143	NA19144	NA19146	NA19147	NA19149	NA19152	NA19153	NA19159	NA19160	NA19171	NA19172	NA19175	NA19184	NA19185	NA19189	NA19190	NA19197	NA19198	NA19200	NA19201	NA19204	NA19206	NA19207	NA19209	NA19210	NA19213	NA19214	NA19222	NA19223	NA19225	NA19235	NA19236	NA19238	NA19239	NA19247	NA19248	NA19256	NA19257	NA19625	NA19700	NA19701	NA19703	NA19704	NA19707	NA19711	NA19712	NA19713	NA19818	NA19819	NA19834	NA19835	NA19900	NA19901	NA19904	NA19908	NA19909	NA19913	NA19914	NA19916	NA19917	NA19920	NA19921	NA19922	NA19923	NA19982	NA19984	NA20126	NA20127	NA20274	NA20276	NA20278	NA20281	NA20282	NA20287	NA20289	NA20291	NA20294	NA20296	NA20298	NA20299	NA20314	NA20317	NA20318	NA20320	NA20321	NA20332	NA20334	NA20339	NA20340	NA20342	NA20346	NA20348	NA20351	NA20355	NA20356	NA20357	NA20359	NA20362	NA20412"
# line_list = line.split()  ### header: NA001 NA002
#
# header_line = '\t'.join(line_list[0:8])
# for col in line_list[9:]:
#     header_line = header_line + '\t' + col + "_1"
#     header_line = header_line + '\t' + col + "_2"
#
# print (header_line)
# manoj_sample_ids_list = []
# for col in sample_ids_list:
#     manoj_sample_ids_list.append(col + "_1")
#     manoj_sample_ids_list.append(col + "_2")
# print(manoj_sample_ids_list)

line = "22	16051249	rs62224609	T	C	100	PASS	.	GT	0|0	0|0	0|0	0|0	0|0	0|0	1|0	0|0	0|0	0|1	0|0	0|1	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|1	0|0	1|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	1|0	0|0	0|0	0|0	1|0	0|0	0|0	0|0	0|0	0|0	1|0	0|0	0|0	0|1	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	1|0	0|0	0|0	1|1	0|0	0|0	0|0	0|0	0|0	0|0	1|0	0|1	0|0	1|0	0|0	1|0	0|0	0|0	0|0	0|0	1|1	1|0	0|0	0|0	1|0	1|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|1	0|0	0|0	0|0	0|0	0|0	0|1	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|1	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|1	0|0	0|0	0|0	0|0	0|0	0|0	0|0	1|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	1|0"
# line = line.split()	###  0|0   0|1
line_list = line.split()

# line_list[9:] = sum([col.split('|') for col in line_list[9:]], [])
# print (line_list[9:])

data_line = '\t'.join(line_list[0:8])
for col in line_list[9:]:
    allele_freq = col.split('|')
    data_line = data_line + '\t' + allele_freq[0]
    data_line = data_line + '\t' + allele_freq[1]
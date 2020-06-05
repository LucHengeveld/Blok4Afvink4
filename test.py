from Bio.Blast import NCBIWWW, NCBIXML

protein = "MGLQARRWASGSRGAAGPRRGVLQLLPLPLPLPLLLLLLLRPGAGRAAAQGEAEAPTLYLWKTGPWGRCMGDECGPGGIQTRAVWCAHVEGWTTLHTNCKQAERPNNQQNCFKVCDWHKELYDWRLGPWNQCQPVISKSLEKPLECIKGEEGIQVREIACIQKDKDIPAEDIICEYFEPKPLLEQACLIPCQQDCIVSEFSAWSECSKTCGSGLQHRTRHVVAPPQFGGSGCPNLTEFQVCQSSPCEAEELRYSLHVGPWSTCSMPHSRQVRQARRRGKNKEREKDRSKGVDPEARELIKKKRNRNRQNRQENKYWDIQIGYQTREVMCINKTGKAADLSFCQQEKLPMTFQSCVITKECQVSEWSEWSPCSKTCHDMVSPAGTRVRTRTIRQFPIGSEKECPEFEEKEPCLSQGDGVVPCATYGWRTTEWTECRVDPLLSQQDKRRGNQTALCGGGIQTREVYCVQANENLLSQLSTHKNKEASKPMDLKLCTGPIPNTTQLCHIPCPTECEVSPWSAWGPCTYENCNDQQGKKGFKLRKRRITNEPTGGSGVTGNCPHLLEAIPCEEPACYDWKAVRLGNCEPDNGKECGPGTQVQEVVCINSDGEEVDRQLCRDAIFPIPVACDAPCPKDCVLSTWSTWSSCSHTCSGKTTEGKQIRARSILAYAGEEGGIRCPNSSALQEVRSCNEHPCTVYHWQTGPWGQCIEDTSVSSFNTTTTWNGEASCSVGMQTRKVICVRVNVGQVGPKKCPESLRPETVRPCLLPCKKDCIVTPYSDWTSCPSSCKEGDSSIRKQSRHRVIIQLPANGGRDCTDPLYEEKACEAPQACQSYRWKTHKWRRCQLVPWSVQQDSPGAQEGCGPGRQARAITCRKQDGGQAGIHECLQYAGPVPALTQACQIPCQDDCQLTSWSKFSSCNGDCGAVRTRKRTLVGKSKKKEKCKNSHLYPLIETQYCPCDKYNAQPVGNWSDCILPEGKVEVLLGMKVQGDIKECGQGYRYQAMACYDQNGRLVETSRCNSHGYIEEACIIPCPSDCKLSEWSNWSRCSKSCGSGVKVRSKWLREKPYNGGRPCPKLDHVNQAQVYEVVPCHSDCNQYLWVTEPWSICKVTFVNMRENCGEGVQTRKVRCMQNTADGPSEHVEDYLCDPEEMPLGSRVCKLPCPEDCVISEWGPWTQCVLPCNQSSFRQRSADPIRQPADEGRSCPNAVEKEPCNLNKNCYHYDYNVTDWSTCQLSEKAVCGNGIKTRMLDCVRSDGKSVDLKYCEALGLEKNWQMNTSCMVECPVNCQLSDWSPWSECSQTCGLTGKMIRRRTVTQPFQGDGRPCPSLMDQSKPCPVKPCYRWQYGQWSPCQVQEAQCGEGTRTRNISCVVSDGSADDFSKVVDEEFCADIELIIDGNKNMVLEESCSQPCPGDCYLKDWSSWSLCQLTCVNGEDLGFGGIQVRSRPVIIQELENQHLCPEQMLETKSCYDGQCYEYKWMASAWKGSSRTVWCQRSDGINVTGGCLVMSQPDADRSCNPPCSQPHSYCSETKTCHCEEGYTEVMSSNSTLEQCTLIPVVVLPTMEDKRGDVKTSRAVHPTQPSSNPAGRGRTWFLQPFGPDGRLKTWVYGVAAGAFVLLIFIVSMIYLACKKPKKPQRRQNNRLKPLTLAYDGDADM"

print("Start Blast...")
result_handle = NCBIWWW.qblast("blastp", "nr", protein, hitlist_size=1)
print("Blasten voltooid.")

with open("XMLTest.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
print("Gegevens in XML bestand gezet")

result_handle = open("XMLTest.xml")

blast_record = NCBIXML.read(result_handle)

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        print(alignment.title)
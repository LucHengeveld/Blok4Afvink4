from flask import Flask, render_template, request
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SearchIO

app = Flask(__name__)


@app.route('/',  methods=["POST", "GET"])
def afvink4():
    if request.method == "POST":

        sequentie = request.form.get("sequentie", "").upper()
        lengte_seq = len(sequentie)

        n_count = sequentie.count("N")
        a_count = sequentie.count("A")
        t_count = sequentie.count("T")
        c_count = sequentie.count("C")
        g_count = sequentie.count("G")
        u_count = sequentie.count("U")

        resultatentext2 = ""
        resultatentext3 = ""

        if lengte_seq == n_count + a_count + t_count + c_count + g_count:
            resultatentext1 = "De ingevoerde sequentie is DNA."

            rna_seq = sequentie.replace("A", "u").replace("T", "a").replace\
                ("C", "g").replace("G", "c").upper()
            resultatentext2 = "De bijbehorende RNA sequentie is: " + rna_seq

            table = {
                'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
                'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
            }
            protein = ""
            if len(sequentie) % 3 == 0:
                for i in range(0, len(sequentie), 3):
                    codon = sequentie[i:i + 3]
                    protein += table[codon]
            resultatentext3 = "De bijbehorende eiwitsequentie is: " + protein

            print("Start Blast...")
            result_handle = NCBIWWW.qblast("blastp", "nr", protein)
            blast_record = NCBIXML.parse(result_handle)
            print(blast_record.hit_def)

        else:
            if lengte_seq == n_count + a_count + u_count + c_count + g_count:
                resultatentext1 = "De ingevoerde sequentie is RNA."

            else:
                if sequentie.isalpha() == True:
                    resultatentext1 = "De ingevoerde sequentie is een eiwit."
                    print("Start Blast...")
                    result_handle = NCBIWWW.qblast("blastp", "nr", sequentie)
                    blast_record = NCBIXML.parse(result_handle)
                    print(blast_record.hit_def)
                    resultatentext3 = "Het meest waarschijnlijke eiwit is: " + blast_record.hit_def

                else:
                    if sequentie.isalpha() == False:
                        resultatentext1 = "Eiwit aan het blasten tegen de " \
                                          "NCBI protein database"
                        return render_template("internetpagina.html",
                                    sequentie=sequentie,
                                    resultatentext1=resultatentext1,
                                    resultatentext2=resultatentext2,
                                    resultatentext3=resultatentext3)

                    resultatentext1 = "De ingevoerde sequentie is geen " \
                                      "DNA, RNA of een Eiwit."

        sequentietxt = "De volgende sequentie is ingevoerd: " + sequentie
        return render_template("internetpagina.html", sequentie=sequentietxt,
                                    resultatentext1=resultatentext1,
                                    resultatentext2=resultatentext2,
                                    resultatentext3=resultatentext3)
    else:
        return render_template("internetpagina.html", sequentie="",
                                    resultatentext1="",
                                    resultatentext2="",
                                    resultatentext3="")


if __name__ == '__main__':
    app.run()

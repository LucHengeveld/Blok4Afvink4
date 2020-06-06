from flask import Flask, render_template, request
from Bio.Blast import NCBIWWW, NCBIXML

app = Flask(__name__)


@app.route('/', methods=["POST", "GET"])
def afvink4():
    """
    :Beschrijving: Geeft de website weer. Haalt de sequentie uit de textbox en
                   berekend hiervan alle gegevens. Nadat deze gegevens zijn
                   berekend wordt dit weergegeven op de website.
    :Return: Webpagina zonder gegevens
             Webpagina met de berekende gegevens (als iets is ingevuld in de
             textbox)
    """
    if request.method == "POST":

        # Gegevens uit textbox halen, dit in caps lock zetten en lengte
        # van de ingevoerde sequentie berekenen.
        sequentie = request.form.get("sequentie", "").upper()
        lengte_seq = len(sequentie)

        # Laat op de website zien welke sequentie is ingevoerd.
        sequentietxt = "De volgende sequentie is ingevoerd: " + sequentie

        # Telt het aantal A, T, C, G, U in de sequentie.
        a_count = sequentie.count("A")
        t_count = sequentie.count("T")
        c_count = sequentie.count("C")
        g_count = sequentie.count("G")
        u_count = sequentie.count("U")

        # Maakt lege strings aan
        resultatentext2 = ""
        resultatentext3 = ""
        resultatentext4 = ""

        # Kijkt op de lengte van de sequentie gelijk is aan A+T+C+G.
        if lengte_seq == a_count + t_count + c_count + g_count:
            resultatentext1 = "De ingevoerde sequentie is DNA."

            # Zet de DNA sequentie om naar RNA om weer te geven op de website.
            rna_seq = sequentie.replace("A", "u").replace("T", "a").replace(
                "C", "g").replace("G", "c").upper()
            resultatentext2 = "De bijbehorende RNA sequentie is: " + rna_seq

            # Kijkt of de ingevoerde sequentie deelbaar is door 3
            # want aminozuren bestaan uit 3 nucleotiden.
            if len(sequentie) % 3 == 0:
                # Berekend de bijbehorende aminozuur sequentie.
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

                # Geeft de aminozuur sequentie door aan de html pagina.
                resultatentext3 = "De bijbehorende aminozuur sequentie is: "\
                                  + protein
            else:
                # Als de sequentie niet door 3 deelbaar is, komt dit op
                # de website te staan
                resultatentext3 = "De ingevoerde DNA sequentie is niet " \
                                  "deelbaar door 3. Hierdoor kunnen niet " \
                                  "alle nucleotiden worden vertaald naar " \
                                  "aminozuren."

            # Blast de RNA sequentie tegen de protein database om het
            # meest waarschijnlijke gen te vinden.
            result_handle = NCBIWWW.qblast("blastx", "nr", rna_seq,
                                           hitlist_size=1)
            # Schrijft de output van BLAST weg in XML bestand.
            with open("XMLTest.xml", "w") as out_handle:
                out_handle.write(result_handle.read())

            # Opent XML bestand en haalt hier de naam uit en geeft het
            # door aan de html pagina.
            result_handle = open("XMLTest.xml", "w")
            blast_record = NCBIXML.read(result_handle)
            for alignment in blast_record.alignments:
                resultatentext4 = "Het meest waarschijnlijke gen is: " + \
                                  alignment.title
                print(resultatentext4)

        else:
            # Kijkt of de sequentie RNA is als het geen DNA is.
            if lengte_seq == a_count + u_count + c_count + g_count:
                resultatentext1 = "De ingevoerde sequentie is RNA."

            else:
                # Kijkt of de sequentie bestaat uit alleen letters als
                # het geen DNA of RNA is.
                if sequentie.isalpha():
                    resultatentext1 = "De ingevoerde sequentie is een eiwit."

                else:
                    # Als het geen DNA, RNA of aminozuur sequentie is,
                    # geeft het dit door aan de website.
                    resultatentext1 = "De ingevoerde sequentie is geen " \
                                      "DNA, RNA of aminozuur sequentie."

        # Return alle gegevens naar de website.
        return render_template("internetpagina.html", sequentie=sequentietxt,
                               resultatentext1=resultatentext1,
                               resultatentext2=resultatentext2,
                               resultatentext3=resultatentext3,
                               resultatentext4=resultatentext4)

    # Return de lege webpagina als nog niets is ingevuld in de textbox
    else:
        return render_template("internetpagina.html", sequentie="",
                               resultatentext1="",
                               resultatentext2="",
                               resultatentext3="",
                               resultatentext4="")


if __name__ == '__main__':
    app.run()

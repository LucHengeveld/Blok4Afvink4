from flask import Flask, render_template, request
import mysql.connector

app = Flask(__name__)


@app.route('/',  methods=["POST", "GET"])
def database():
    if request.method == "POST":
        resultatenlijst = []
        zoekopdracht = request.form.get("zoek", "")

        conn = mysql.connector.connect(host="ensembldb.ensembl.org",
                                       user="anonymous",
                                       db="homo_sapiens_core_95_38")
        cursor = conn.cursor()
        query = "select * from gene where description like \"%" + zoekopdracht\
               + "%\""
        cursor.execute(query)

        for rij in cursor:
            resultatenlijst.append(rij[9])
        cursor.close()
        conn.close()

        rangeresultatentext = "Alleen de eerste 20 resultaten worden " \
                              "weergegeven. Dit kan worden aangepast " \
                              "onderaan de body in de range bij de for loop."

        return render_template("internetpagina.html",
                               resultaten=resultatenlijst,
                               resultatentext=
                               'Gevonden resultaten:',
                               resultatentextrange=rangeresultatentext)
    else:
        return render_template('internetpagina.html',
                               resultaten='',
                               resultatentext='',
                               resultatentextrange='')


if __name__ == '__main__':
    app.run()

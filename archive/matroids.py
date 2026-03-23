arango = pyimport("arango")
client = ArangoClient(hosts="https://matroid.mathematik.uni-siegen.de:443")
db = client.db("MatroidsDB", username="matroid", password="matroid")
matroids = db.collection("matroids_split_public")


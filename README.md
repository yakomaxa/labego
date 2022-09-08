# labego
show ABEGO by label on PyMOL (label + abego = labego)

# Usage
In PyMOL,

```run labego.py```

```labego (target)```

labego only shows ABEG part of ABEGO

```labegO (target)```

labegO also shows O of ABEGO, but this requires full backbone atoms (works fine for computational models such as Rosetta and AlphaFold-derived models)

```lapsego (target)```

lapsego and lapsegO divides B-region into two sub-regions named P and B, so it's called laPSego(O)

See this [paper](https://www.jstage.jst.go.jp/article/biophysico/18/0/18_bppb-v18.017/_html/-char/en)
 if you wonder why such sub-division is helpful

# License
MIT License

Author Koya Q Sakuma 2020

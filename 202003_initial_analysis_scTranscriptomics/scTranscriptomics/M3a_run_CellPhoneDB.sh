source ~/cpdb-venv/bin/activate

cellphonedb database generate --result-path database --user-gene ~/farm/CellPhoneDB-data_smallmolecules/data/gene_input_all.csv --user-complex ~/farm/CellPhoneDB-data_smallmolecules/data/sources/complex_curated.csv --user-interactions ~/farm/CellPhoneDB-data_smallmolecules/data/sources/interaction_curated.csv

cellphonedb method analysis ./meta.tsv ./counts.csv --database ./out/database/cellphonedb_user_2019-11-29-15_57.db --counts-data hgnc_symbol --output-path ./out/ --threshold 0.1
# conda activate ersilia

import os
from ersilia import ErsiliaModel

RESULTS_FOLDER = "../../results/"

model_ids = ["eos2r5a", "eos7pw8"]

data_folder = "../../data/generated/zairachem"

for input_file in os.listdir(data_folder):
    if not input_file.endswith(".csv"): continue
    for model_id in model_ids:
        output_file = os.path.abspath(os.path.join(RESULTS_FOLDER, "ersilia-{0}-{1}".format(model_id, os.path.split(input_file)[-1])))
        input_file = os.path.abspath(os.path.join(data_folder, input_file))
        if os.path.exists(output_file):
            continue
        print(input_file, output_file)
        em = ErsiliaModel(model=model_id)
        em.serve()
        em.api(input=input_file, output=output_file)
        em.close()
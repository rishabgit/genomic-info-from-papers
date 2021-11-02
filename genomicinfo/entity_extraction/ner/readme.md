# How to train NER:
Create baseline model trained on [IDP4+](https://www.tagtog.net/jmcejuela/IDP4plus/-settings) corpus by following these steps:  
  
1.

```   
>> from genomicinfo.entity_extraction.ner.data_prep import NERDataPrep 
>> NERDataPrep('IDP4+')
```

2. 

python train_ner.py --model_name_or_path dmis-lab/biobert-base-cased-v1.1 --train_file training_data/tagtog_IDP4+_anndoc/train_dev.json --validation_file training_data/tagtog_IDP4+_anndoc/devel.json --text_column_name tokens --label_column_name tags --pad_to_max_length --max_length 192 --per_device_train_batch_size 8 --learning_rate 2e-5 --num_train_epochs 10 --output_dir saved_models/nala --seed 1


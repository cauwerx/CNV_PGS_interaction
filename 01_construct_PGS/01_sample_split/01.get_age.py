# Obtain year of birth information 

################################################################################
# Prologue apache spark
################################################################################

# Load libraries
import pyspark 
import dxpy
import dxdata
import pandas as pd

# Spark initialization
sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)

# Get dispensed dataset ID
dispensed_dataset_id = dxpy.find_one_data_object(
    typename = "Dataset",
    name = "app*.dataset",
    folder = "/",
    name_mode = "glob"
)["id"]

dataset = dxdata.load_dataset(id=dispensed_dataset_id)

################################################################################
# STEP 1: Retrieve participant information
################################################################################

# Access main participant entity
participant = dataset["participant"]

# p34 is year of birth
field_names = ["eid", "p34"]

# Extracting fields
spark_df = participant.retrieve_fields(names=field_names,
                                       coding_values="replace",
                                       engine=dxdata.connect())

# Convert to pandas
df = spark_df.toPandas()
# Check I have no NaN values
assert not df.isna().any().any(), "NaN values found!"
# eid and date of birth should be int
df = df.astype(int)

# Rename columns
df.rename(columns={'p34': 'birth_year'}, inplace=True)
df.rename(columns={'eid': 'IID'}, inplace=True)

################################################################################
# STEP 2: Save participant information
################################################################################

df.to_parquet('age.all.parquet', engine='pyarrow') 

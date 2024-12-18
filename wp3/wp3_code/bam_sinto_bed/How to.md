# BAM File Processing and Fragment File Creation

This folder contains the code to **add barcodes** to the BAM files and **create fragment files** from the extended BAM files.

The subfolder **`twoVMs`** contains two scripts to process the **NAPKON BAM files** on two virtual machines (VMs):

- **B.** uses the script: `bam_sinto_bed_vm1.sh`  
- **S.** uses the script: `bam_sinto_bed_vm2.sh`  


## Steps to Process the BAM Files

### 0. Obtain the Environment File  
Download the `environment.yml` file from the Git repository and save it on the VM.

---

### 1. Create the Conda Environment  
Run the following command to create the Conda environment:  

```bash
conda env create -f environment.yml
```
### 2. Activate the Environment
Activate the newly created environment:
```bash
conda activate wp3
```
### 3. Execute the Job Script
Run the corresponding job script for your VM.
```bash
nohup bash ./wp3/wp3_code/bam_sinto_bed/twoVMs/bam_sinto_bed_vm2.sh > output2.log 2>&1 &
```

### Troubleshooting
In case the code breaks, most errors will be written into the error files located in the error directory.


## Merging Fragment Files
Once the FRAGMENTS files are created, they will be merged by B. into one consolidated file named:

*napkon_fragments.bed*



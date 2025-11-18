# Instructions for Globus Transfer

# 1. Install Globus
**You only need to do this once**. If you already have set up everything from a previous run, please go to section #3.

a. First, obtain a [globus id](https://www.globusid.org/create).<br>
b. Navigate to a directory on the server and install globus connect personal with `wget`. I navigate to my scratch directory (`/lustre/isaac24/scratch/netid`), but you can go anywhere you'd like.
```
wget https://downloads.globus.org/globus-connect-personal/linux/stable/globusconnectpersonal-latest.tgz
``` 
c. Extract the files from the downloaded tarball.
```
tar xzf globusconnectpersonal-latest.tgz # this will produce a versioned globusconnectpersonal directory
cd globusconnectpersonal-x.y.z # replace x.y.z with the version number you see
```
d. Start up the program to initialize.
```
./globusconnectpersonal
```
This will interactively guide you through the set-up process. It produces a link; copy this into your browser of choice, copy the resulting code, and paste into the terminal. It will also ask you to name your local endpoint; choose whatever name you like. Set up is now complete!

# 2. Set up local endpoint
a. Start the program.
```
./globusconnectpersonal -start &
```
b. Download the globus CLI. Use conda or mamba, as it is better at resolving dependencies (i.e., less errors!)
```
conda create -n globus
conda activate globus
conda install -c conda-forge globus-cli
```
c. Use the new conda install of globus to identify your local endpoint.
```
globus endpoint local-id
```
When I run this command, I got a long string of characters and numbers: 76aeec42-9348-11f0-8647-0e840c2393b5 (yours will be different). To make it easier, you can also save this ID as a shortcut, which will make the transfer command a lot easier. Below, I named my local endpoint 'allyson':
```
export allyson="76aeec42-9348-11f0-8647-0e840c2393b5"
```
# 3. Updating permissions for a local endpoint
When using Globus, you will only be able to transfer files to and from directories that are set to accessible. To add appropriate directories, make an edit to the config file:
```
nano ~/.globusonline/lta/config-paths
```
If never run before, you will see something like this:
```
~/,0,1
```
* The first entry is the absolute path of desired directory (in this case, it's the home directory).<br>
* The second entry is the sharing flag. `0` disables sharing and `1` enables sharing.<br>
* The third entry is the R/W flag, which enables write-access. `0` allows read-only access and `1` allows write-access.<br>

**Be sure to add every directory you want to transfer files to**. You can repeat this step as many times as needed, especially if your target directory changes between projects. Once a directory is added, it will remain saved. Make sure all required directories are added before starting the transfer to avoid errors.

# 4. Initiating the transfer
a. Find the endpoint where your files are stored.
```
globus endpoint search 'NAME'
```
For example, let's say I want to find the directory where Makhali's Smokies files are. Usually the UT Genomics Core will include the name of the PI. So, I could search the following:
```
globus endpoint search 'Owings'
```
This search gives me the following list:
```
ID                                   | Owner                                                        | Display Name                                   
------------------------------------ | ------------------------------------------------------------ | -----------------------------------------------
67dc70ba-18ca-11ed-a909-fd3165076336 | owe5924@northwestern.edu                                     | Turash_lab_PC                                  
56bf1d31-4709-4116-b1f4-a44b4fb9ef94 | oleary@brown.edu                                             | BrownU_OwenLeary_GDrive                        
e1d3f7bf-b05a-48c5-b3d8-d4c21dba584d | subed042@umn.edu                                             | DPR data from Owen                             
0088e0dd-3108-4357-8bee-6137fd5de42d | 4636e469-e85a-470e-a923-6bea93357f6b@clients.auth.globus.org | Imperial College London - David Owen - NGS Data
30699fd8-9992-46bf-a43d-3f09e39d0a35 | jbemmels@computecanada.ca                                    | Kelp ARGs Owens Lab                            
d369cfeb-d64d-4c01-ba8b-52486031c704 | sbsbioinformatics@globusid.org                               | NGS Data - Imperial College London - David Owen
9504e1f4-72ae-11ea-9610-0afc9e7dd773 | 4c984b40-a0b2-4d9e-b132-b32735905e23@clients.auth.globus.org | nigelmic_owenfunk                              
780c0a2a-1ce5-457c-8837-d65b0f65c2da | 4636e469-e85a-470e-a923-6bea93357f6b@clients.auth.globus.org | Nucleome Therapeutics - NGS Data Lauren Owens  
81df261f-52fc-40f8-949e-2a03dfe4e751 | statius9@ufl.edu                                             | oweiss_Nicholas_data_dir                       
5c54d3a0-2c82-4a02-b5e3-312bdcb81ee1 | subed042@umn.edu                                             | Owens data                                     
1139aa12-8734-411f-8463-f1ec24a90ce3 | subed042@umn.edu                                             | Owens data                                     
0ba67447-ece0-40c3-abcd-2a18090b00f0 | jmill165@tennessee.edu                                       | Owings UTK Illumina Data 20251107              
3499691d-44b9-4139-b952-eddf956e3e24 | rmhzmh0@identity.6a4e18.8540.data.globus.org                 | rd0216_nick_owen                               
1d18e9ca-a156-4148-b161-0ff566d14896 | viannam@ufl.edu                                              | SC_genomes_Owen                                
47e55769-be0a-4a2e-8d5d-1fe81ceacf5c | vlarivie@computecanada.ca                                    | CO.SHS
```
The directory I want would be `Owings UTK Illumina Data 20251107` and the corresponding ID is `0ba67447-ece0-40c3-abcd-2a18090b00f0`. Let's save that as a shortcut to make life easier.
```
export owings_illumina_nov2025="0ba67447-ece0-40c3-abcd-2a18090b00f0"
```
Now we can initiate the transfer. Here is the skeletal code to initiate a transfer:
```
globus transfer initial_endpoint:/path/to/files my_endpoint:/desired/path/to/files --recursive --label "Custom Label"
```
Let's give an actionable example. Makhali's files are in a directory called `fastq` in the `Owings UTK Illumina Data 20251107` directory. Let's say I want to transfer these files to my personal scratch directory (`/lustre/isaac24/scratch/adekovic`) using my local endpoint (nickname `allyson` or 76aeec42-9348-11f0-8647-0e840c2393b5).

Here would be an example command:
```
globus transfer owings_illumina_nov2025:/fastq allyson:/lustre/isaac24/scratch/adekovic --recursive --label "Smokies Metagenomics Samples"
```
* `owings_illumina_nov2025` and `allyson` are the shortcuts for the long endpoint IDs. Alternatively, you can paste the entire ID in this section and it will work just fine.<br>
* `--recursive` makes sure everything within the `fastq` directory gets copied.<br>
* `--label` creates a custom label for the job; this will show up in the web version/any email notifications about progress.

Now you have your samples! Happy coding :)

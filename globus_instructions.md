# Instructions for Globus Transfer

# 1. Install Globus and set up local endpoint.
**You only need to do this once**. If you already have set up everything from a previous run, please go to section #2.

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

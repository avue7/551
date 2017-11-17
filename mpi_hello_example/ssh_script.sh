#!/bin/bash

echo -e "Now running the ssh-agent script for MPI jobs."
ssh-agent bash
ssh-add ~/.ssh/id_rsa
cat ssh_script
echo -e "Now will ask you fo your passphrase."


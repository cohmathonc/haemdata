#!/bin/bash
watch 'squeue --Format=jobid:8,username:12,name:70,state:8,timeused:10,timeleft:10,starttime,endtime,numcpus:5,minmemory:8,partition:8,nodelist:13,account:8 -u $USER --sort -N'


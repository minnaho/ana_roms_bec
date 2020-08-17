#!/bin/csh -f
#  roms.cmd
#
#  UGE job for roms built Mon Aug 17 10:40:05 PDT 2020
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/home/m/minnaho/ana_roms_bec/roms.joblog.$JOB_ID
#$ -o /u/home/m/minnaho/ana_roms_bec/roms.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/home/m/minnaho/ana_roms_bec/roms
#  arguments       = roms.in
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Resources requested
#
#$ -l h_data=4000M,h_rt=1:00:00
# #
#  Name of application for log
#$ -v QQAPP=job
#  Email address to notify
#$ -M minnaho@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
#
# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "job serial"
  set qqidir    = /u/home/m/minnaho/ana_roms_bec
  set qqjob     = roms
  set qqodir    = /u/home/m/minnaho/ana_roms_bec
  cd     /u/home/m/minnaho/ana_roms_bec
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for roms built Mon Aug 17 10:40:05 PDT 2020"
  echo ""
  echo "  roms directory:"
  echo "    "/u/home/m/minnaho/ana_roms_bec
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
#
  echo ""
  echo "roms started on:   "` hostname -s `
  echo "roms started at:   "` date `
  echo ""
#
# Run the user program
#
  source /u/local/Modules/default/init/modules.csh
  module load intel/13.cs
#
  echo roms "roms.in" \>\& roms.output.$JOB_ID
  echo ""
  /usr/bin/time /u/home/m/minnaho/ana_roms_bec/roms roms.in >& /u/home/m/minnaho/ana_roms_bec/roms.output.$JOB_ID
#
  echo ""
  echo "roms finished at:  "` date `
#
# Cleanup after serial execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/home/m/minnaho/ana_roms_bec/roms.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/job.log.serial
  if (`wc -l /u/home/m/minnaho/ana_roms_bec/roms.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
	head -50 /u/home/m/minnaho/ana_roms_bec/roms.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
	echo " "  >> /u/local/apps/queue.logs/job.log.serial
	tail -10 /u/home/m/minnaho/ana_roms_bec/roms.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  else
	cat /u/home/m/minnaho/ana_roms_bec/roms.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  endif
  exit (0)

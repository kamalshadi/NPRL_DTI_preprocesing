'''
Preprocessing pipeline for NPRL DTI dataset written based on
HCP minimal preprocessing guidelines for diffusion MRIs

usage:
 -> preprocess -i <study_folder> -t <fsl_readout_time>

Requirement:
Python 2.7

Url:
https://github.com/kamalshadi/NPRL_DTI_preprocesing

Author:
Kamal Shadi

Contact:
kshadi3@gatech.edu

Date:
2017/06/17

Revision:
1.0


Copyright
Neural Plasticity Research Laboratory (NPRL) - Emory University
'''

import os
import argparse
import sys
import copy
import numpy as num
from math import cos,sin
from shutil import copyfile
import json



#################### Setting input arguments ###########################
parser = argparse.ArgumentParser(description=\
'Pipeline to preprocess DTI data')
parser.add_argument('-i', action='store', dest='root_folder',\
help='Root folder containing a subfolder called "Inputs" having'+\
' all subjects\' raw images.',\
required=True)
parser.add_argument('-t', action='store', dest='fsl_readout',\
help='Total readout time (FSL definition). you can get this from'+\
' MRIconvert tool.',\
type = float, default=0.0508)



############ decompressing nifti for exploreDTI compatability ##########
def gzip(fd):
	print fd
	for (dp, dns, fns) in os.walk(fd):
		for w in fns:
			if '.nii' in w and '.gz' in w:
				fn = dp + w
				os.system("gunzip -d "+fn)


########################### bvec file utilities ########################
def read_bvec(fn):
	with open(fn) as f:
		st = f.readlines()
	k = {}
	for i,w in enumerate(st):
		k[i] = w.split()
	return k

def read_bval(fn):
	with open(fn) as f:
		st=f.read()
	return [xx.strip() for xx in st.split()]


def rot_bvec(cur,theta_x,theta_y,theta_z):
	Rx = num.zeros((3,3))
	Ry = num.zeros((3,3))
	Rz = num.zeros((3,3))
	Rx[0,0] = cos(theta_x)
	Rx[0,1] = sin(theta_x)
	Rx[1,0] = -sin(theta_x)
	Rx[1,1] = cos(theta_x)
	Rx[2,2] = 1.0

	Ry[0,0] = cos(theta_y)
	Ry[0,1] = sin(theta_y)
	Ry[1,0] = -sin(theta_y)
	Ry[1,1] = cos(theta_y)
	Ry[2,2] = 1.0

	Rz[0,0] = cos(theta_z)
	Rz[0,1] = sin(theta_z)
	Rz[1,0] = -sin(theta_z)
	Rz[1,1] = cos(theta_z)
	Rz[2,2] = 1.0

	R = num.dot(num.dot(Rx,Ry),Rz)
	_R = num.linalg.inv(R)
	return num.dot(_R,cur)


####################### Cleaning input directory #######################
def clean_copy(fd):
	if fd[-1]!='/':
		fd = fd + '/'
	fn = fd + 'data_eddy.eddy_parameters'
	with open(fn) as f:
		st = f.readlines()

	bvec = read_bvec(fd+'data.bvec')
	out_bvec = copy.copy(bvec)
	cur = num.zeros((3,1))
	for q,w in enumerate(st):
		w = w.split()
		tmp = [float(xx) for xx in w]
		theta_x = tmp[3]
		theta_y = tmp[4]
		theta_z = tmp[5]
		cur[0] = bvec[0][q]
		cur[1] = bvec[1][q]
		cur[2] = bvec[2][q]
		out = rot_bvec(cur,theta_x,theta_y,theta_z)
		out_bvec[0][q] = out[0][0]
		out_bvec[1][q] = out[1][0]
		out_bvec[2][q] = out[2][0]
	with open(fd+'data_rotated.bvec','w') as f:
		for i in range(3):
			cur = out_bvec[i]
			st = ' '.join([str(xx) for xx in cur])
			st = st + '\n'
			f.write(st)



################### Main pipeline script ###############################
def pipeline(fd,fsl_readout):
	out = fd.split('/Inputs/')[0]+'/Outputs'
	if not os.path.isdir(out):
		os.mkdir(out)
	sub_out = fd.replace('/Inputs/','/Outputs/')
	if not os.path.isdir(sub_out):
		os.mkdir(sub_out)
	#~ print ';;;;'
	if fd[-1]!='/':
		fd = fd + '/'
	merging = True
	ac = True
	topup = True
	avg = True
	bet = True
	bvec = True
	bval = True
	index = True
	eddy = True
	clr = True

	#finding files
	b0_up = None
	b0_dn = None
	raw_up = None
	raw_dn = None
	bvec_up = None
	bvec_dn =None
	for (dp, dns, fns) in os.walk(fd):
		for w in fns:
			if b0_rl[0] in w and b0_rl[1] in w and '.nii' in w:
				b0_up = dp + w
			elif b0_lr[0] in w and b0_lr[1] in w and '.nii' in w:
				b0_dn = dp + w
			elif dwi_rl[0] in w and dwi_rl[1] in w and '.nii' in w:
				raw_up = dp + w
			elif dwi_lr[0] in w and dwi_lr[1] in w and '.nii' in w:
				raw_dn = dp + w
			elif dwi_rl[0] in w and dwi_rl[1] in w and '.bvec' in w:
				bvec_up = dp + w
			elif dwi_lr[0] in w and dwi_lr[1] in w and '.bvec' in w:
				bvec_dn = dp + w
			elif dwi_rl[0] in w and dwi_rl[1] in w and '.bval' in w:
				bval_up = dp + w
			elif dwi_lr[0] in w and dwi_lr[1] in w and '.bval' in w:
				bval_dn = dp + w
			else:
				pass
		break

	if b0_up is None or b0_dn is None or raw_up is None\
	or raw_dn is None \
	or bvec_up is None or bvec_dn is None:
		print 'Error in input files of directory '+dp
		print 'Make sure you provide all needed files for all subjects.'
		sys.exit(1)

	if merging:
		# forming nodif
		print 'Juxtaposing the data in the right order...'
		with open(bval_up) as f:
			tmp = f.read()
			tmp = tmp.strip()
		with open(bval_up,'w') as f:
			f.write(tmp)
		with open(bval_dn) as f:
			tmp = f.read()
			tmp = tmp.strip()
		with open(bval_dn,'w') as f:
			f.write(tmp)

		os.system('select_dwi_vols '+raw_up+ ' '+bval_up+' '+fd+'nodif_up1 0')
		os.system('select_dwi_vols '+raw_dn+ ' '+bval_dn+' '+fd+'nodif_dn1 0')
		os.system('select_dwi_vols '+raw_up+ ' '+bval_up+' '+fd+'weighted_up 1000')
		os.system('select_dwi_vols '+raw_dn+ ' '+bval_dn+' '+fd+'weighted_dn 1000')
		os.system('fslmerge -t '+fd+'nodif_up '+b0_up+' ' \
		+fd+'nodif_up1')
		os.system('fslmerge -t '+fd+'nodif_dn '+b0_dn+' '\
		+fd+'nodif_dn1')
		os.system('fslmerge -t '+fd+'nodif '+fd+'nodif_up '+fd+\
		'nodif_dn')

		# forming data volumes

		os.system('fslmerge -t '+fd+'data '+fd+'weighted_up '+\
		fd+'weighted_dn')
		os.system('fslmerge -t '+fd+'data '+fd+'nodif '+fd+'data')
	if ac:
		c1 = C['volumes']['b0_rl']+C['volumes']['dwi_rl'][0]
		c2 = C['volumes']['b0_lr']+C['volumes']['dwi_lr'][0]
		# forming acqparams
		print 'Making acqparams...'
		with open(fd+'acqparams.txt','w') as f:
			for i in range(c1):
				f.write('1 0 0 '+str(fsl_readout)+'\n')
			for i in range(c2):
				f.write('-1 0 0 '+str(fsl_readout)+'\n')
	if topup:
		# running topup
		print 'Running topup...'
		cmd = '''topup --imain=IMAIN --config=b02b0.cnf --out=OUT \
		--iout=IOUT --datain=DATAIN'''
		cmd = cmd.replace('IMAIN',fd+'nodif')
		cmd = cmd.replace('IMAIN',fd+'data')
		cmd = cmd.replace('IOUT',fd+'b0_topup_brain')
		cmd = cmd.replace('OUT',fd+'b0_topup')
		cmd = cmd.replace('DATAIN',fd+'acqparams.txt')
		os.system(cmd)

	# average nodif brain

	if avg:
		print 'Making average brain...'
		cmd = 'fslmaths '+fd+'b0_topup_brain -Tmean '+fd+'_brain'
		os.system(cmd)



	# extracting brain

	if bet:
		print 'Making brain mask...'
		cmd = 'bet '+fd+'_brain '+fd+'brain -m -f 0.3'
		os.system(cmd)


	# forming bvec and bval
	if bvec:
		print 'Forming bvec...'
		up = read_bvec(bvec_up)
		dn = read_bvec(bvec_dn)
		bup = read_bval(bval_up)
		bdn = read_bval(bval_dn)
		with open(fd+'data.bvec','w') as f:
			for i in range(3):
				tmp = ['0.0']*(C['volumes']['b0_rl']+C['volumes']['dwi_rl'][0]+\
				C['volumes']['b0_lr']+C['volumes']['dwi_lr'][0])
				for q,item in enumerate(up[i]):
					if int(bup[q])==0:
						continue
					else:
						tmp = tmp + [item]
				for q,item in enumerate(dn[i]):
					if int(bdn[q])==0:
						continue
					else:
						tmp = tmp + [item]
				st = ' '.join(tmp)
				f.write(st+'\n')




	# forming bval
	if bval:
		print 'Forming bval...'
		tmp = ['0']*(C['volumes']['b0_rl']+C['volumes']['dwi_rl'][0]+\
		C['volumes']['b0_lr']+C['volumes']['dwi_lr'][0])
		tmp = tmp + ['1000']*(C['volumes']['dwi_rl'][1]+\
		C['volumes']['dwi_lr'][1])
		st = ' '.join(tmp)
		with open(fd+'data.bval','w') as f:
			f.write(st)


	# forming index
	if index:
		print 'Forming index...'
		n0_up=C['volumes']['b0_rl']+C['volumes']['dwi_rl'][0]
		n0_dn=C['volumes']['b0_lr']+C['volumes']['dwi_lr'][0]
		n_up=C['volumes']['dwi_rl'][1]
		n_dn=C['volumes']['dwi_lr'][1]
		tmp = range(1,n0_up+n0_dn+1)
		tmp = [str(xx) for xx in tmp]
		tmp = tmp + [str(n0_up)]*n_up + [str(n0_up+n0_dn)]*n_dn
		tmp = ' '.join(tmp)
		with open(fd+'index','w') as f:
			f.write(tmp)

	# running eddy
	if eddy:
		print 'running eddy...'
	#eddy --index=index.txt --bvecs=data.bvecs
	#--bvals=data.bval --mask=brain
	#--acqp=acqparam.txt --topup=b0_topup --out=data_eddy
		cmd = '''eddy --imain=IMAIN --index=INDEX --bvecs=BVEC
		 --bvals=BVAL --mask=MASK --topup=TOPUP --out=OUT --acqp=ACQ'''
		cmd = cmd.replace('INDEX',fd+'index')
		cmd = cmd.replace('BVAL',fd+'data.bval')
		cmd = cmd.replace('BVEC',fd+'data.bvec')
		cmd = cmd.replace('MASK',fd+'brain_mask')
		cmd = cmd.replace('ACQ',fd+'acqparams.txt')
		cmd = cmd.replace('TOPUP',fd+'b0_topup')
		cmd = cmd.replace('OUT',fd+'data_eddy')
		cmd = cmd.replace('IMAIN',fd+'data')
		os.system(cmd)
	#print fd
	clean_copy(fd)
	print 'Decompressing the .gz to .nii files...'
	gzip(fd)
	print 'Copying results to output folder...'
	files = ['data_eddy.nii','nodif.nii','brain.nii','brain_mask.nii',\
	'data.bval','data_rotated.bvec']
	for cur in files:
		copyfile(fd+cur,sub_out+'/'+cur)
	return 1





	# saving results to output
	#corrected data
	#bval, bvec, acqparams, index
	if clr:
		pass

results = parser.parse_args()


################### Input settings from config.json #############################
rt = results.root_folder
if rt[-1] != '/':
	rt = rt + '/'
with open(rt+'config.json') as f:
    C = json.load(f)
try:
	b0_rl = C['filenames']['b0_rl']
except KeyError:
	b0_rl = None
try:
	b0_lr = C['filenames']['b0_lr']
except KeyError:
	b0_lr = None

try:
	dwi_rl = C['filenames']['dwi_rl']
except KeyError:
	print 'Error in config file:Please include dwi_rl key for right to left (or up down) dwi images'
	sys.exit(1)


try:
	dwi_lr = C['filenames']['dwi_lr']
except KeyError:
	print 'Error in config file:Please include dwi_lr key for left to right (or down up) dwi images'
	sys.exit(1)

####################### Starting the pipeline ##########################
dirc = rt+'Inputs/'
if not os.path.isdir(dirc):
	print 'Error:directory '+dirc+' not found!'
	sys.exit()
for (dp, dns, fns) in os.walk(dirc):
	if dp[0] == '/':
		dp = '/'+dp.strip('/')+'/'
	else:
		dp = dp.strip('/')+'/'
	for i,w in enumerate(dns):
		sub = dp+w
		print '------------------------------------------------'
		print str(i+1)+'/'+str(len(dns))
		print 'Processing starts for subject: '+ sub
		pipeline(sub,results.fsl_readout)

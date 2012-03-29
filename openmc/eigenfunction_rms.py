#!/usr/bin/python
# Filename: eigenfunction_rms.py

import h5py
import numpy as np
import cPickle
import matplotlib.pyplot as plt
import sys
#
class EigenFunction:
#
 	'''Represents the reference eigenfunction'''
 	reference = 0.
#
	def __init__(self,data):
#
		'''Initializes the eigenfunction'''
		self.function = 0.
		self.data = data
		self.meanfunction = 0.
#
        def read_hdf5(self,h5_file,cycle):
#
		'''Read data from HDF5 file'''
                self.cycle = cycle
                f = h5py.File(h5_file,'r')
                group = '/cycle'+str(cycle)+'/'+self.data
                dataset = f[group]
                self.function = np.empty(dataset.shape,dataset.dtype)
                dataset.read_direct(self.function)
		self.function = self.function
                self.iamref = 'F'
#
	def set_reference(self):
#
		'''Sets instance to be reference calc'''
		self.iamref = 'T'
		EigenFunction.reference = self.function
#
	def compute_rms(self):
#
		'''Computes RMS value'''
		Np = self.function.size
		Np = 177 
                tmp = (self.meanfunction - EigenFunction.reference)**2
                tmp2 = tmp.sum()
                self.rms = np.sqrt((1.0/float(Np))*tmp2)
#
	def compute_all_rms(self,function):
#
		'''Computes RMS value'''
		Np = self.function.size
		Np = 177
		tmp = (function - EigenFunction.reference)**2
		tmp2 = tmp.sum()
		rms = np.sqrt((1.0/float(Np))*tmp2)
		return rms
#
def read_runs(runpath,hdfile,cycle_start,cycle_end,run_start,run_end,data):
	runlist = []
	tmp = EigenFunction(data)
	tmp.read_hdf5(runpath+str(run_start)+'/'+hdfile,cycle_start)
	indices = tmp.function.shape # extent of all dimensions
	i = cycle_start
	while i <= cycle_end:
		meantmp = EigenFunction(data) # init mean object
		j = run_start
		runs = np.zeros((run_end - run_start + 1,indices[0],indices[1],indices[2],indices[3])) # init runs array
		while j <= run_end:
			tmp.read_hdf5(runpath+str(j)+'/'+hdfile,i) # read hdf5 file
			runs[j-run_start] = tmp.function # put function into runs
			j += 1
		meantmp.meanfunction = np.average(runs, axis=0) # compute the mean
		meantmp.function = runs
		runlist.append(meantmp)
		print 'Read in from path: '+runpath+' Cycle: '+str(i)
		i += 1
	return runlist

def create_reference(runpath,hdfile,cycle,run_start,run_end,data):

#	 calculate reference solution
	print 'Calculating Reference solution...'
	tmp = EigenFunction(data)
	tmp.read_hdf5(runpath+str(run_start)+'/'+hdfile,cycle) # load first eigenfunction
	print 'Read in: '+runpath+str(1)+hdfile
	indices = tmp.function.shape # extent of all dimensions
	ref = np.zeros((run_end,indices[0],indices[1],indices[2],indices[3])) # initialize ref array
	ref[0] = tmp.function # set the first run in ref
	i = run_start + 1
	while i <= run_end: # begin loop around all runs
		tmp.read_hdf5(runpath+str(i)+'/'+hdfile,cycle)
		print 'Read in: '+runpath+str(i)+hdfile
		ref[i - run_start] = tmp.function
		i += 1
	meanref = np.average(ref, axis=0) # compute average of all runs
	EigenFunction.reference = meanref # set to global space in EigenFunction instances
	return meanref

def compute_rms(runlist):

	# calculate rms array
	print 'Calculating rms...'
	rms = np.zeros(len(runlist))
	i = 0
	while i < len(runlist):
		runlist[i].compute_rms()
		rms[i] = runlist[i].rms
		i += 1
	return rms

def plot_rms(rms_1mil,rms_cmfd,rms_64mil):
	print 'Generating plot...'
	ax = plt.subplot(111)
	size_1 = rms_1mil.shape[0]
	x_1 = np.linspace(1,size_1,size_1)*1e6
	y_1 = (rms_1mil[0]/x_1[0]**(-0.5))*x_1**(-0.5)
	size_64 = rms_64mil.shape[0]
	x_64 = np.linspace(1,size_64,size_64)*64e6
	plt.loglog(x_1,rms_1mil*100,'b+')
	plt.loglog(x_1,y_1*100,'k-')
        plt.loglog(x_1,rms_cmfd*100,'g--')
	plt.loglog(x_64,rms_64mil*100,'ro')
	ax.xaxis.grid(True,'minor')
	ax.yaxis.grid(True,'minor')
	ax.xaxis.grid(True,'major',linewidth=2)
	ax.yaxis.grid(True,'major',linewidth=2)
	plt.xlabel('# of Total Neutron Histories (active cycles)')
	plt.ylabel('RMS Error [%]')
	plt.legend(('1 million (10 runs)','Ideal Error 1 mil','CMFD 1mil','64 million (10 runs)'))
	return

def plot_source(source):
	plt.figure()
	X = np.linspace(1,17,17)
	Y = np.linspace(1,17,17)
	Y,X = np.meshgrid(Y,X)
	plt.contourf(X,Y,source[0,:,:,0],50)
	plt.colorbar()
	plt.xlabel('Mesh Cell in x-direction')
	plt.ylabel('Mesh Cell in y-direction')
	plt.title('OPR Converged Fission Source Distribution')
#
def plot_all_rms(onemil,cmfd,sixfourmil):

	# new plot
	plt.figure()
	ax = plt.subplot(111)
	ax.xaxis.grid(True,'minor')
	ax.yaxis.grid(True,'minor')
	ax.xaxis.grid(True,'major',linewidth=2)
	ax.yaxis.grid(True,'major',linewidth=2)
	
	# loop and plot onemil
	size = onemil[0].function.shape[0]
	rmsmat_1 = np.zeros([size,len(onemil)])
	i = 0
	while i < size:
		j = 0
		while j < len(onemil):
			rmsmat_1[i,j]=onemil[j].compute_all_rms(onemil[j].function[i])
			j += 1
		x = np.linspace(1,len(onemil),len(onemil))*1e6
		plt.loglog(x,rmsmat_1[i]*100,'b+')
		print 'Computing Run '+str(i)+'...'
		i += 1
        # loop and plot onemil
        size = cmfd[0].function.shape[0]
        rmsmat_cmfd = np.zeros([size,len(cmfd)])
        i = 0
        while i < size:
                j = 0
                while j < len(cmfd):
                        rmsmat_cmfd[i,j]=cmfd[j].compute_all_rms(cmfd[j].function[i])
                        j += 1
                x = np.linspace(1,len(cmfd),len(cmfd))*1e6
                plt.loglog(x,rmsmat_cmfd[i]*100,'g--')
                print 'Computing Run '+str(i)+'...'
                i += 1

	# loop and plot sixfourmil
	size = sixfourmil[0].function.shape[0]
	rmsmat_64 = np.zeros([size,len(sixfourmil)])
	i = 0
	while i < size:
		j = 0
		while j < len(sixfourmil):
			rmsmat_64[i,j] = sixfourmil[j].compute_all_rms(sixfourmil[j].function[i])
			j += 1
		x = np.linspace(1,len(sixfourmil),len(sixfourmil))*64e6
		plt.loglog(x,rmsmat_64[i]*100,'ro')
		i += 1

	# show all plots 
        plt.xlabel('# of Total Neutron Histories (active cycles)')
        plt.ylabel('RMS Error [%]')
        plt.show()
#
if __name__ == "__main__":

	if sys.argv[1] == 'restart':

		# load in data
		print 'Loading input...'
		filein = open('rms.out','r')
		output = cPickle.load(filein)
		filein.close()
                rms_onemil = output['1milrms']
                rms_sixfourmil = output['64milrms']
		meanref = output['ref']
		EigenFunction.reference = meanref

		# plot rms
	        plot_rms(rms_onemil,rms_sixfourmil)

		# plot mean source distribution
		plot_source(meanref)

	elif sys.argv[1] == 'interactive':

                # load in data
		print 'Loading input...'
		filein = open('rms.out','r')
		output = cPickle.load(filein)
		filein.close()
		rms_onemil = output['1milrms']
		rms_sixfourmil = output['64milrms']
		meanref = output['ref']
		EigenFunction.reference = meanref

		# pop an interactive python shell
	 	from IPython import embed
		embed()

	else:

		# calculate reference solution
		runpath = '/media/Backup/opr_runs/64mil/run'
		hdfile = 'output.h5'
		cycle = 210
		run_start = 1
		run_end = 25 
		data = 'openmc_src'
		meanref = create_reference(runpath,hdfile,cycle,run_start,run_end,data)

		# calculate rms for 1 million case
                runpath = '/media/Backup/opr_runs/1mil/run'
                hdfile = 'output.h5'
		cycle_start = 201
		cycle_end = 840
		run_start = 1
		run_end = 10 
		data = 'openmc_src'
		onemil = read_runs(runpath,hdfile,cycle_start,cycle_end,run_start,run_end,data)
		rms_onemil = compute_rms(onemil)

                # calculate rms for 1 million case
                runpath = '/media/Backup/opr_runs/1mil/run'
                hdfile = 'output.h5'
                cycle_start = 201
                cycle_end = 840
                run_start = 1
                run_end = 10
                data = 'cmfd_source'
                cmfd = read_runs(runpath,hdfile,cycle_start,cycle_end,run_start,run_end,data)
                rms_cmfd = compute_rms(cmfd)

		# calculate rms for 64 million case
		runpath = '/media/Backup/opr_runs/64mil/run'
		hdfile = 'output.h5'
		cycle_start = 201
		cycle_end = 210
		run_start = 1
		run_end = 10 
		data = 'openmc_src'
		sixfourmil = read_runs(runpath,hdfile,cycle_start,cycle_end,run_start,run_end,data)
                rms_sixfourmil = compute_rms(sixfourmil)

		# write out numpy array to binary file
		print 'Writing output...'
		output = {}
		output.update({'1milrms':rms_onemil})
		output.update({'64milrms':rms_sixfourmil})
		output.update({'ref':meanref})
		fileout = open('rms.out','wb')
		cPickle.dump(output,fileout)
		fileout.close()

		# plot rms
		plot_rms(rms_onemil,rms_cmfd,rms_sixfourmil)

		# plot mean source distribution
		plot_source(meanref)

		# plot all rms
		plot_all_rms(onemil,cmfd,sixfourmil)

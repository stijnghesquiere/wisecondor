##############################################################################
#                                                                            #
#    Plot test results of WISECONDOR.                                        #
#    Copyright(C) 2013  TU Delft & VU University Medical Center Amsterdam    #
#    Author: Roy Straver, r.straver@vumc.nl                                  #
#                                                                            #
#    This file is part of WISECONDOR.                                        #
#                                                                            #
#    WISECONDOR is free software: you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    WISECONDOR is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with WISECONDOR.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                            #
##############################################################################


import sys
import pickle
import gcc

import numpy
import math
import cutoff

import matplotlib
import argparse

#    Output adapted by Stijn Ghesquiere, s.ghesquiere@mumc.nl, stijn@applesnail.net
def plotResults(sample, markedBins, kept, kept2, outputFile, zScoresDict,zSmoothDict,blindsDict,sampleName):
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib.collections import BrokenBarHCollection
	from matplotlib.backends.backend_pdf import PdfPages # stijn addition for multipage pdf output

	# stijn: find largest chromosome in terms of number of bins
	def largest_chromosome():
		bins_largest_chrom = 0
		for chrom in range(1,23):
			if int(len(sample[str(chrom)])) > bins_largest_chrom:
				bins_largest_chrom = int(len(sample[str(chrom)])) 
		#print 'largest is:', bins_largest_chrom
		return(bins_largest_chrom)

	colorSample			= 'blue'
	colorReference 		= 'red'
	colorMarkedBin		= (0.5,1,0.35)
	colorMarkedSure		= (0.3,0.6,0.07)#'green'
	colorMarkedSmooth	= (1,0.6,1)
	colorMarkedSureSm	= (0.4,0.2,0.4)
	colorHorzHelper		= (0.7,0.7,0.7)
	colorHorzMarker		= 'orange'
	colorBlinds			= (0.85,0.85,0.85)
	colorWaste			= (0.7,0.7,0.7)

	binScalar = 320.0 # = 2x DPI 160 from document --> PDF
	edgeSize = 0.15
	grid_columns = largest_chromosome()

	wastedBins = dict()
	for chrom in range(1,23):
		wastedBins[str(chrom)] = []
		for bin in range(len(sample[str(chrom)])-1):
			wastedBins[str(chrom)].append(len(getReference(sample,str(chrom),bin,[],0,4,1)) <= 3)

	def ideogram(chrom):
		color_lookup = {
			'gneg': (1., 1., 1.),
			'gpos25': (.6, .6, .6),
			'gpos50': (.4, .4, .4),
			'gpos75': (.2, .2, .2),
			'gpos100': (0., 0., 0.),
			'acen': (.8, .4, .4),
			'gvar': (.8, .8, .8),
			'stalk': (.9, .9, .9),
			'select': (0., 0., 1.)
		}
		fin = open('./cytoBand.txt')
		prevChrom = ''
		for line in fin:
			lchrom, start, stop, label, stain = line.strip().split('\t')
			if ('chr'+ str(chrom)) == lchrom:
				start = int(start)/float(binSize)
				stop = int(stop)/float(binSize)				
				width = stop - start
				center = start + width/2.
				plt.axvline(x=center,  ymin=0-edgeSize, ymax=1, linewidth=width*(binScalar/int(len(sample[str(chrom)]))*2), color=color_lookup[stain])
				plt.text(center, 0.2, label, fontsize=6,backgroundcolor='white',verticalalignment='center',horizontalalignment='center',rotation='vertical')
			elif prevChrom != lchrom and prevChrom == chrom:
				break
			prevChrom = lchrom
		fin.close()
		frame1 = plt.gca()
		frame1.axes.get_xaxis().set_visible(False)
		frame1.axes.get_yaxis().set_visible(False)

	def printMarkers(chrom,grid_columns):
		#binWidth = 1/float(len(sample[str(chrom)]))*binScalar 
		binWidth =  binScalar/grid_columns * 2.0
		zSmooth = zSmoothDict[str(chrom)]
		for bin in range(len(zSmooth)):
			if abs(zSmooth[bin]) > 3:
				plt.axvline(x=bin, ymin=1-edgeSize*2, ymax=1, linewidth=binWidth, color=colorMarkedSmooth)

		for region in kept2:
			if str(region[0]) == chrom:
				for bin in range(region[1],region[2]+1):
					plt.axvline(x=bin, ymin=1-edgeSize, ymax=1, linewidth=binWidth, color=colorMarkedSureSm)

		marks = [1] * len(sample[str(chrom)])
		for bin in markedBins:
			if str(bin[0]) == chrom:
				plt.axvline(x=bin[1], ymin=0, ymax=edgeSize*2, linewidth=binWidth, color=colorMarkedBin)

		for region in kept:
			if str(region[0]) == chrom:
				for bin in range(region[1],region[2]+1):
					plt.axvline(x=bin, ymin=0, ymax=edgeSize, linewidth=binWidth, color=colorMarkedSure)

	def drawBlinds(chrom,grid_columns):
		binWidth =  binScalar/grid_columns * 2
		for bin in blindsDict[str(chrom)]:
			plt.axvline(x=bin, linewidth=binWidth, color=colorBlinds)
		for bin in range(len(wastedBins[str(chrom)])):
			if wastedBins[str(chrom)][bin]:
				plt.axvline(x=bin, linewidth=binWidth, color=colorWaste)

	def preparePlot(chrom,grid_columns):
		drawBlinds(chrom,grid_columns)
		printMarkers(str(chrom),grid_columns)
	print 'Plotting Z-Scores'
	plt.figure(2)




	def overviewPlot():
		fig = plt.figure(figsize=(11.7, 8.3),dpi=160) # Stijn: new page in pdf
		matplotlib.use('PDF')
		for chrom in range(1,23):
			plt.subplot2grid( (69,grid_columns), (chrom*3-3,0), rowspan=2, colspan = int(len(sample[str(chrom)])) )
			preparePlot(chrom,grid_columns)
			zSmooth = zSmoothDict[str(chrom)]
			zKeep = []
			for val in zScoresDict[str(chrom)]:
				if not math.isnan(val):
					zKeep.append(val)
			zTotal = numpy.sum(zKeep) / numpy.sqrt(len(zKeep))

			plt.axhline(y=0, linewidth=0.75, color=colorHorzHelper) # stijn: z-score helper lines
			plt.axhline(y=3, linewidth=0.50, color=colorHorzHelper)
			plt.axhline(y=-3, linewidth=0.50, color=colorHorzHelper)

			plt.plot(zSmooth,colorReference)
			plt.plot(zScoresDict[str(chrom)],colorSample)

			plt.ylabel(chrom)
			plt.xlim(0,len(sample[str(chrom)])-2)
			plt.ylim(-12,12)
			frame1 = plt.gca()
			frame1.axes.get_xaxis().set_visible(False)

			for tick in frame1.axes.get_yaxis().get_major_ticks():
				tick.label.set_fontsize(3)
		pdf_pages.savefig(fig) # Stijn: close page 


	def PlotChromosome(chrom,y_scale,grid_columns):
		preparePlot(chrom,grid_columns)
		zSmooth = zSmoothDict[str(chrom)]
		zKeep = []
		for val in zScoresDict[str(chrom)]:
			if not math.isnan(val):
				zKeep.append(val)
		zTotal = numpy.sum(zKeep) / numpy.sqrt(len(zKeep))
		if y_scale != 0:
			plt.ylim(y_scale*(-1),y_scale) # Stijn: fixed Y axis scale
		plt.axhline(y=0, linewidth=0.75, color=colorHorzHelper) # stijn: z-score helper lines
		plt.axhline(y=3, linewidth=0.50, color=colorHorzHelper)
		plt.axhline(y=-3, linewidth=0.50, color=colorHorzHelper)
		plt.plot(zSmooth,colorReference)
		plt.plot(zScoresDict[str(chrom)],colorSample)
		plt.plot(sample[str(int(chrom))],'black')
		
		plt.ylabel(chrom)
		plt.xlim(0,len(sample[str(chrom)])-2)

		frame1 = plt.gca()
		frame1.axes.get_xaxis().set_visible(False)


	pdf_pages = PdfPages(outputFile+'.zscore.pdf') # multipage pdf

	overviewPlot()

	for chrom in [13,18,21]:
		fig = plt.figure(figsize=(11.7, 8.3),dpi=160) # Stijn: new page in pdf
		plt.subplot2grid( (4,int(len(sample[str(chrom)]))), (0,0), rowspan=1, colspan = int(len(sample[str(chrom)])) )
		print '------>', int(len(sample[str(chrom)])) 
		PlotChromosome(chrom,6,int(len(sample[str(chrom)])))

		# plot with fixed Y-axis (z-score)
		plt.subplot2grid( (4,int(len(sample[str(chrom)]))), (1,0), rowspan=1, colspan = int(len(sample[str(chrom)])) )
		PlotChromosome(chrom,12,int(len(sample[str(chrom)])))

		# plot with variable Y-axis (z-score)
		plt.subplot2grid( (4,int(len(sample[str(chrom)]))), (2,0), rowspan=1, colspan = int(len(sample[str(chrom)])) )
		PlotChromosome(chrom,0,int(len(sample[str(chrom)])))

		# ideogram - chromosomes
		plt.subplot2grid( (4,int(len(sample[str(chrom)]))), (3,0), colspan = int(len(sample[str(chrom)])) )
		ideogram(chrom)

		pdf_pages.savefig(fig) # Stijn: close page

	pdf_pages.close() # close pdf



if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Plot results generated by WISECONDOR',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('plotdata', type=str,
					   help='datafile as output by test.py')
	parser.add_argument('outfile', type=str,
					   help='output file to store plot in, .pdf is added as extension')
					   
	parser.add_argument('-mpluse', default='agg', type=str, 
					help='make matplotlib use another backend for plotting')

	args = parser.parse_args()

	print '# Script information:'
	print '\n# Settings used:'
	matplotlib.use(args.mpluse)
	argsDict = args.__dict__
	argsKeys = argsDict.keys()
	argsKeys.sort()
	for arg in argsKeys:
		print '\t'.join([arg,str(argsDict[arg])])
		
	print '\n# Processing:'
	print 'Loading:\tSample:\t' + args.plotdata
	sampleData = pickle.load(open(args.plotdata,'rb'))
	
	sampleName=args.plotdata.split("/")[-1].split(".")[0]
			
	plotResults( \
		sampleData['sample'], \
		sampleData['markedBins'], \
		sampleData['kept'], \
		sampleData['kept2'], \
		args.outfile, \
		sampleData['zScoresDict'], \
		sampleData['zSmoothDict'], \
		sampleData['blindsDict'], \
		sampleData['wastedBins'], \
		sampleName
		)
	print '\n# Finished'

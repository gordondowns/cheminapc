# -*- coding: utf-8 -*-
"""
Created on Wed Mar 07 09:23:08 2018

@author: gordi
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import minimize
from time import time
#from CifFile import ReadCif

#import warnings
#warnings.filterwarnings("error")

class XrayPhase(object):
    """The x-ray diffraction data of one material. Can be either a Profile (full pattern XY data)
    or a Dif (peaks only, which will be broadened). Has scalings, RIR, chemistry, and name of corresponding file."""
    
    def __init__(self,file_path=None):
        self.RIR = None
        self.scaling_guess = None
        self.scaling_optimized = None
        self.scaling_bounds = (0.0,None)
        self.refine_scaling = True
        self.chemistry = None
        self.relative_contribution = None
        self.absolute_contribution = None
        self.file_path = file_path
        if file_path is None:
            self.file_nickname = None
        else:
            self.file_nickname = file_path.split('/')[-1].split('\\')[-1].split('.')[0]
    
    def SetScalingGuess(self,scaling_guess,set_both=True):
        self.scaling_guess = scaling_guess
        if set_both:
            self.scaling_optimized = scaling_guess
    
    def SetScalingOptimized(self,scaling_optimized):
        self.scaling_optimized = scaling_optimized

    def SetScalingBounds(self,scaling_bounds):
        self.scaling_bounds = scaling_bounds
        
    def FitScaling(self,boolean):
        self.refine_scaling = boolean

    def GetArea(self,start_2theta=None,end_2theta=None):
        pass
    
    def GetAreaSloppy(self,start_2theta,end_2theta):
        pass
    
    def ReadFile(self,file_path): 
        pass

    def SetRelativeContribution(self,percentage):
        self.relative_contribution = percentage
    def GetRelativeContribution(self):
        return self.relative_contribution
    def SetAbsoluteContribution(self,percentage):
        self.absolute_contribution = percentage
    def GetAbsoluteContribution(self):
        return self.absolute_contribution



class XrayProfile(XrayPhase):
    """The x-ray diffraction data of one material, with full pattern XY data. Can be read in from MDI, XY, or CSV."""
    
    def __init__(self, file_path, is_input_profile=False, take_every_nth_point=1, twotheta_ranges=[(0.0,90.0)], twotheta_offset=0.0, print_warnings=True):
        super(XrayProfile,self).__init__(file_path)
        self.xy_data = None
        self.xy_data_for_plotting = None
        self.x_data_interpolated = None
        self.y_data_interpolated = None
        self.area_with_scaling_1 = None
        if(is_input_profile):
            self.SetScalingGuess(1.0)
        self.ReadFile(file_path,is_input_profile,take_every_nth_point,twotheta_ranges,twotheta_offset,print_warnings=print_warnings)
    
    def GetInterpolatedXYData(self, xValues):
        #print "GetInterpolatedXYData.xValues",xValues
        
        self.x_data_interpolated = np.array(xValues)
        self.y_data_interpolated = np.interp(xValues,self.xy_data[0],self.xy_data[1])
        #FIXME: implement this
        return np.array([self.x_data_interpolated,self.y_data_interpolated])
    
    def GetScaledPattern(self,xValues=None,scaling=None,true_for_optimized_false_for_guess=True):
        if(scaling is None):
            if true_for_optimized_false_for_guess:
                scaling = self.scaling_optimized
            else:
                scaling = self.scaling_guess
        #if(not xValues is None): #FIXME: implement this
        self.GetInterpolatedXYData(xValues)
        return np.multiply( self.y_data_interpolated, scaling )
    
    def GetArea(self, xValues= None, scaling=None, start_2theta=None, end_2theta=None):
        from scipy import integrate
        if scaling is None:
            scaling = self.scaling_optimized
        if xValues is not None:
            self.GetInterpolatedXYData(xValues)
            self.area_with_scaling_1 = integrate.simps(self.y_data_interpolated, xValues)
        
        return self.area_with_scaling_1 * scaling
        #FIXME: Implement ability to set start and end 2thetas to None
    
    
    def ReadFile(self, file_path, is_input_profile=False, take_every_nth_point=1, twotheta_ranges=[(0.0,90.0)], twotheta_offset=0.0, print_warnings=True):
        
        RIR_found, chemistry_found = False, False
        RIR = 0.0
        chemistry = {}
        should_read_next_line = True
        mdi_header_found, xy_data_found = False, False
        
        
        f = open(file_path, 'r')
#        print f.readlines()
        while(True):
            if should_read_next_line:
                line = f.readline()
            else:
                should_read_next_line = True
#            print '-'
#            print line
            if line == "":
                break
            elif "".join(line.split()).startswith("#RIR"):
                RIR_found = True
                RIR = float(line.lstrip().split()[2])
            elif "".join(line.split()).startswith("#CHEMISTRY"):
                chemistry_found = True
                while(True):
                    line = f.readline()
                    should_read_next_line = False
                    linesplit = line.split()
                    if len(linesplit) != 2 or not linesplit[1].replace('.','').isdigit():
                        # chemistry not found on this line
                        break
                    else:
                        # chemistry found on this line
                        name,percentage = linesplit
                        chemistry[name] = float(percentage)
            else:
                try: # look for MDI format
                    linesplit = line.replace('.','').split()
#                    if line.startswith('2.0 0.05'):
#                        pass
#                        print linesplit
                    if (len(linesplit) == 7 and linesplit[0].isdigit() and linesplit[1].isdigit() and
                    linesplit[5].isdigit() and linesplit[6].isdigit()):
                        
                        mdi_header_found = True
                        linesplit = line.split()
                        x_start = float(linesplit[0])
                        #x_increment = float(line2split[1])
                        #wavelength = float(linesplit[4])
                        x_end = float(linesplit[5])
                        x_num = int(linesplit[6])
                        
                        x_data = np.linspace(x_start,x_end,num=x_num)
                        y_data = []
                        
                        while(True):
                            line = f.readline()
                            should_read_next_line = False
                            linesplit = line.replace('.','').split()
                            if "".join(linesplit).isdigit():
                                y_data += [float(i) for i in line.split()]
                            else:
                                break
                except:
                    pass
                try: # look for XY or CSV format
                    linesplit = line.replace(',',' ').replace('.','').replace('e','').replace('E','').replace('-','').replace('+','').split()
                    
                    if len(linesplit) == 2 and linesplit[0].isdigit() and linesplit[1].isdigit():
                        linesplit = line.replace(',',' ').split()
                        x_data = [float(linesplit[0])]
                        y_data = [float(linesplit[1])]
                        xy_data_found = True
                        
                        while(True):
                            line = f.readline()
                            should_read_next_line = False
                            linesplit = line.replace(',',' ').replace('.','').replace('e','').replace('E','').replace('-','').replace('+','').split()
                            if len(linesplit) == 2 and linesplit[0].isdigit() and linesplit[1].isdigit():
                                linesplit = line.replace(',',' ').split()
                                x_data += [float(linesplit[0])]
                                y_data += [float(linesplit[1])]
                            else:
                                break
                except:
                    pass            
        if not (mdi_header_found and len(x_data) == len(y_data) or xy_data_found):
            f.close()
            raise SyntaxError("No data found in profile: '"+file_path+"'\n"
                             +'Profile files must contain either XY data or MDI data.\n'
                             +'- XY files must contain at least one line that is\n'
                             +'  a pair of numbers separated by whitespace and/or a comma.\n'
                             +'- MDI files must contain standard MDI header line of the form:\n'
                             +'  x_start x_increment scan_rate xray_source xray_wavelength x_end number_of_points')
            
        f.close()
        # At this point, at least x_data and y_data have been found, and RIR and chemistry have been looked for.
        # Done parsing the file. Time to process the data.
        
        self.xy_data_unprocessed = np.array([x_data,y_data])
        
        # apply 2theta offset: add offset to every 2theta value
        if twotheta_offset != 0.0:
            x_data = [x + twotheta_offset for x in x_data]
        
        # set minimum intensity to 0.0 (if not input profile) and set maximum intensity to 100.0
        maxY = max(y_data)
        if is_input_profile:
            y_data = np.array([y_data[i]*100.0/maxY for i in range(len(y_data))])
        else:
            minY = min(y_data)
            y_data = np.array([(y_data[i]-minY) *100.0/(maxY-minY) for i in range(len(y_data))])
        
        self.xy_data_for_plotting = np.array([x_data,y_data])

        if take_every_nth_point > 1:
            x_data = [x_data[i] for i in range(len(x_data)) if i % take_every_nth_point == 0]
            y_data = [y_data[i] for i in range(len(y_data)) if i % take_every_nth_point == 0]
        
        # remove data points outside of desired range
        if twotheta_ranges != [(0.0,90.0)]:
            new_x_data,new_y_data = [],[]
            for i in range(len(x_data)):
                for pair in twotheta_ranges:
                    if x_data[i] >= pair[0] and x_data[i] <= pair[1]:
                        new_x_data.append(x_data[i])
                        new_y_data.append(y_data[i])
                        break # It is sufficient for the data to be in only one interval. So, go to next x_data.
            x_data,y_data = new_x_data,new_y_data
        
        self.xy_data = np.array([x_data,y_data])
        self.RIR = RIR
        self.chemistry = chemistry
        if print_warnings:
            if(not is_input_profile and not RIR_found):
                print("Warning: no RIR found in amorphous file: "+file_path)
            if(not chemistry_found):
                print("Warning: no chemistry found in amorphous file: "+file_path)
        
    def StripPeaks(self, keep = "amorphous", num_loops=300, delta_bounds=(3,7), graph = False):
        import random
#        from scipy.interpolate import UnivariateSpline
        
        background_part = np.array(self.xy_data[1])
        
        for j in range(num_loops):
            delta = random.randint(delta_bounds[0],delta_bounds[1])
            temp_background_part = np.array(background_part)
            for i in range(delta,len(background_part)-delta):
                background_part[i] = min(temp_background_part[i], (temp_background_part[i-delta]+temp_background_part[i+delta])/2)
        
#        s = UnivariateSpline(self.xy_data[0],background_part,s=1)
#        background_part = s(self.xy_data[0])
        
        
        crystalline_part = self.xy_data[1] - background_part
                
        if graph:
            plt.figure(1,figsize=(30, 15))
            plt.subplot(222)
            plt.plot(self.xy_data[0],self.xy_data[1],label="raw data")
            plt.plot(self.xy_data[0],background_part,label="background/amorphous component")
            plt.plot(self.xy_data[0],crystalline_part,label="peaks/crystalline component")
            plt.title("peak stripping for profile "+self.file_path)
            plt.grid(True)
            plt.legend(bbox_to_anchor=(1.0, -0.05))
            plt.xlabel('2-theta (deg)')
            plt.ylabel('intensity')
            plt.axis([min(self.xy_data[0]),max(self.xy_data[0]), 0.0 ,20.0])
            plt.show()
        
        if keep == "amorphous" or keep == "background":
            self.xy_data = np.array([self.xy_data[0],background_part])
        elif keep == "crystalline" or keep == "peaks":
            self.xy_data = np.array([self.xy_data[0],crystalline_part])
        
        return crystalline_part, background_part
    
class XrayDif(XrayPhase):
    """The x-ray diffraction data of one material, with only peak XY data.
    Reads only difs from AMCSD or RRUFF (no file modifications needed)."""
    
    def __init__(self, file_path=None, print_warnings=True):
        super(XrayDif,self).__init__(file_path)
        self.FWHM_guess = None
        self.FWHM_optimized = None
        self.FWHM_bounds = (None,None)
        self.refine_FWHM = True
        self.pseudovoigt_parameter_guess = None
        self.pseudovoigt_parameter_optimized = None
        self.pseudovoigt_parameter_bounds = (0.0,1.0)
        self.refine_pseudovoigt_parameter = False
        self.peak_data = None
        self.HKLs = None
        self.d_spacings = None
        self.cell_parameters_original = None #for copypasta
        self.cell_parameters_guess = None
        self.cell_parameters_optimized = None
        self.cell_parameter_bounds = [(0.1,None),(0.1,None),(0.1,None),(30.0,None),(30.0,None),(30.0,None)]
        self.refine_cell_parameters = True # equivalent to [True,True,True,True,True,True] # haven't implemented this yet
        self.x_ray_wavelength = None
        self.space_group = None
        self.num_unique_cell_parameters = None
        self.crystal_system = None
        self.ReadFile(file_path,print_warnings)
    
    def SetFWHMGuess(self,FWHM_guess,set_both=True):
        self.FWHM_guess = FWHM_guess
        if set_both:
            self.FWHM_optimized = FWHM_guess
    
    def SetFWHMOptimized(self,FWHM_optimized):
        self.FWHM_optimized = FWHM_optimized
    
    def SetFWHMBounds(self,FWHM_bounds):
        self.FWHM_bounds = FWHM_bounds
        
    def FitFWHM(self,boolean):
        self.refine_FWHM = boolean
    
    def SetPseudoVoigtParameterGuess(self,pseudovoigt_parameter_guess,set_both=True):
        self.pseudovoigt_parameter_guess = pseudovoigt_parameter_guess
        if set_both:
            self.pseudovoigt_parameter_optimized = pseudovoigt_parameter_guess
    
    def SetPseudoVoigtParameterOptimized(self,pseudovoigt_parameter_optimized):
        self.pseudovoigt_parameter_optimized = pseudovoigt_parameter_optimized
        
    def SetPseudoVoigtParameterBounds(self,pseudovoigt_parameter_bounds):
        self.pseudovoigt_parameter_bounds = pseudovoigt_parameter_bounds
        
    def FitPseudoVoigtParameter(self,boolean):
        self.refine_pseudovoigt_parameter = boolean
    
    def SetCellParametersGuess(self, cell_parameters_guess):
        self.cell_parameters_guess = cell_parameters_guess
    
    def SetCellParametersOptimized(self, cell_parameters_optimized):
        self.cell_parameters_optimized = cell_parameters_optimized
    
    def SetCellParameterBoundsUsingValues(self, cell_parameter_bounds):
        self.cell_parameter_bounds = cell_parameter_bounds
    
    def SetCellParameterBoundsUsingPercentage(self, percentage = 5.0, true_for_optimized_false_for_guess = True):
        if true_for_optimized_false_for_guess:
            self.cell_parameter_bounds = [(self.cell_parameters_optimized[i]*(1.0-percentage/100.0),self.cell_parameters_optimized[i]*(1.0+percentage/100.0)) for i in range(6)]
        else:
            self.cell_parameter_bounds = [(self.cell_parameters_guess[i]*(1.0-percentage/100.0),self.cell_parameters_guess[i]*(1.0+percentage/100.0)) for i in range(6)]

    def SetCellParameterBoundsUsingDifferences(self, plusminus = [0.05,0.05,0.05,1.0,1.0,1.0], true_for_optimized_false_for_guess = True):
        if true_for_optimized_false_for_guess:
            self.cell_parameter_bounds = [(self.cell_parameters_optimized[i]-plusminus[i],self.cell_parameters_optimized[i]+plusminus[i]) for i in range(6)]
        else:
            self.cell_parameter_bounds = [(self.cell_parameters_guess[i]-plusminus[i],self.cell_parameters_guess[i]+plusminus[i]) for i in range(6)]

    def FitCellParameters(self, boolean):
        self.refine_cell_parameters = boolean
    
    def GetArea(self, FWHM=None, scaling=None, start_2theta=None, end_2theta=None):
        total_area = 0.0
        
        if FWHM is None:
            FWHM = self.FWHM_optimized
        if scaling is None:
            scaling = self.scaling_optimized
        
        sigma = FWHM / 2.35482
        
        if (start_2theta is not None) and (end_2theta is not None):
            for i in range(len(self.peak_data[0])):
                total_area += self.peak_data[1][i] * (norm.cdf(x=end_2theta,loc=self.peak_data[0][i],scale=sigma) - norm.cdf(x=start_2theta,loc=self.peak_data[0][i],scale=sigma))
        
        if (start_2theta is not None) and (end_2theta is None):
            for i in range(len(self.peak_data[0])):
                total_area += self.peak_data[1][i] * (1.0 - norm.cdf(x=start_2theta,loc=self.peak_data[0][i],scale=sigma))
        
        if (start_2theta is None) and (end_2theta is not None):
            for i in range(len(self.peak_data[0])):
                total_area += self.peak_data[1][i] * norm.cdf(x=end_2theta,loc=self.peak_data[0][i],scale=sigma)
        
        if (start_2theta is None) and (end_2theta is None):
            for i in range(len(self.peak_data[0])):
                total_area += self.peak_data[1][i] # area under a Gaussian is 1
        
        return total_area * scaling
    

    
    def GetBroadenedPattern(self,xValues,FWHM=None,scaling=None,pseudovoigt_parameter=None,true_for_optimized_false_for_guess=True):
        #turns the peaks into full xy data for plotting
        
        eta = pseudovoigt_parameter
        xin = self.peak_data[0]
        yin = self.peak_data[1]
        
        if (FWHM is None):
            if true_for_optimized_false_for_guess:
                FWHM = self.FWHM_optimized
            else:
                FWHM = self.FWHM_guess
        if (scaling is None):
            if true_for_optimized_false_for_guess:
                scaling = self.scaling_optimized
            else:
                scaling = self.scaling_guess
        if (eta is None):
            if true_for_optimized_false_for_guess:
                eta = self.pseudovoigt_parameter_optimized
            else:
                eta = self.pseudovoigt_parameter_guess
        elif eta > 1.0 or eta < 0.0:
            print("Pseudo-voigt parameter eta must be between 0 and 1.")
            print("The given value of eta is "+str(eta))
            return


        k = 0.01
        
        # convenient Gaussian parameter is sigma
        sigma = FWHM / 2.35482
        tempG = 1/(sigma * np.sqrt(2 * np.pi))
        
        yout = np.zeros_like(xValues)
        
        if eta == 0.0: #then pure Gaussian
        
            for j in range(0,len(xin)): # jth peak of this dif
                
                if(yin[j] < 0.001):
                    continue
                
                delta = np.sqrt( abs( -2.0*sigma**2 *np.log( max(k/(yin[j]*scaling*tempG), 0.0000001) ) ) )
                
                i_start,i_end = np.searchsorted(xValues,[xin[j]-delta,xin[j]+delta])
                
                for i in range(i_start,i_end): # ith x value
                    yout[i] += yin[j] * scaling * tempG * np.exp( - (xValues[i] - xin[j])**2 / (2 * sigma**2))
        
        else: # then pseudo-voigt
            # convenient Lorentzian parameter is gamma
            gamma = FWHM / 2.0
            tempL = gamma / np.pi
            
            for j in range(0,len(xin)): # jth peak of this dif
            
                if(yin[j] < 0.001):
                    continue
                
                # Gaussian component
                delta = np.sqrt( abs( -2.0*sigma**2 *np.log( max(k/(yin[j]*scaling*tempG), 0.0000001) ) ) )
                # Lorentzian component
                delta += np.sqrt( abs( -2.0*sigma**2 *np.log( max(k/(yin[j]*scaling*tempG), 0.0000001) ) ) )
                # FIXME ^
                for i in range(len(yout)): # ith x value
#                    try:
                        # Gaussian component
                        yout[i] += yin[j] * scaling * (1.0-eta) * tempG * np.exp( - (xValues[i] - xin[j])**2 / (2 * sigma**2))
                        # Lorentzian component
                        yout[i] += yin[j] * scaling * eta * tempL / ( (xValues[i]-xin[j])**2+gamma**2 )
#                    except RuntimeWarning:
#                        print 'sigma', sigma
#                        print 'scaling',scaling
#                        print 'tempG',tempG
#                        print 'tempL', tempL
#                        print 'eta',eta
#                        #AAAAAAAAAAAAAAAAAA
                        
#                for i in range(len(yout)): # ith x value
#                    # Gaussian component
#                    yout[i] += yin[j] * scaling * (1.0-eta) * tempG * np.exp( - (xValues[i] - xin[j])**2 / (2 * sigma**2))
#                    # Lorentzian component
#                    yout[i] += yin[j] * scaling * eta * tempL / ( (xValues[i]-xin[j])**2+gamma**2 )

#        else:
#            h = scaling
##            x0 = xin[j]
#            f = FWHM
#            lFraction = eta
#            gFraction = 1.0-lFraction
#            M_LN2 = np.log(2.0)
#            
#            gamma = FWHM / 2.0
#            sigma = sigma
#
##            // Lorentzian parameter gamma...fwhm/2
##            g = f / 2.0;
##            double gSquared = g * g;
##                        
##            // Gaussian parameter sigma...fwhm/(2*sqrt(2*ln(2)))...gamma/sqrt(2*ln(2))
##            double sSquared = gSquared / (2.0 * M_LN2);
#            gSquared = (f/2.0)**2
#            sSquared = gSquared/(2.0*M_LN2)
#            
#            gSquared = gamma**2
#            sSquared = sigma**2
#            
##            for (size_t i = 0; i < nData; ++i) {
##                double xDiffSquared = (xValues[i] - x0) * (xValues[i] - x0);
##            
##                out[i] = h * (gFraction * exp(-0.5 * xDiffSquared / sSquared) +
##                   (lFraction * gSquared / (xDiffSquared + gSquared)));
##               }
#            
#            yout = np.zeros_like(xValues)
#            for j in range(len(xin)):
#                yout += np.array([yin[j] * (gFraction * np.exp(-0.5 * (x-xin[j])**2 / sSquared) + (lFraction * gSquared / ((x-xin[j])**2 + gSquared)))for x in xValues])
            

        return np.array(yout)

    
    
    def GetPseudoVoigtPattern(self,xValues,FWHM=None,scaling=None,pseudovoigt_parameter=None,true_for_optimized_false_for_guess=True):
        #turns the peaks into full xy data as pseudo-voigts rather than Gaussians
        
        # eta is the parameter that determines how much Gaussian and how much Lorentzian there is.
        # eta = 0 means completely Gaussian and eta = 1 is completely Lorentzian
        # V(x,FWHM) = eta*L(x,FWHM) + (1-eta)*G(x,FWHM) with 0 < eta < 1
        
        xin = self.peak_data[0]
        yin = self.peak_data[1]
        
        eta = pseudovoigt_parameter
        
        if (FWHM is None):
            if true_for_optimized_false_for_guess:
                FWHM = self.FWHM_optimized
            else:
                FWHM = self.FWHM_guess
        if (scaling is None):
            if true_for_optimized_false_for_guess:
                scaling = self.scaling_optimized
            else:
                scaling = self.scaling_guess
        if (eta is None):
            if true_for_optimized_false_for_guess:
                scaling = self.pseudovoigt_parameter_optimized
            else:
                scaling = self.pseudovoigt_parameter_guess
        elif eta > 1.0 or eta < 0.0:
            print("Pseudo-voigt parameter eta must be between 0 and 1.")
            print("The given value of eta is "+str(eta))
            return
 
        k = 0.01
        
        # convenient Gaussian parameter is sigma
        sigma = FWHM / 2.35482
        # convenient Lorentzian parameter is gamma
        gamma = FWHM / 2
        
        tempG = 1/(sigma * np.sqrt(2 * np.pi))
        tempL = gamma / np.pi
        
        yout = np.zeros_like(xValues)
        
        for j in range(0,len(xin)): # jth peak of this dif
            
            if(yin[j] < 0.001):
                continue
            
            for i in range(len(xin)): # ith x value
                # Gaussian component
                yout[i] += yin[j] * scaling * (1.0-eta) * tempG * np.exp( - (xValues[i] - xin[j])**2 / (2 * sigma**2))
                # Lorentzian component
                yout[i] += yin[j] * scaling * eta * tempL / ( (xValues[i]-xin[j])**2+gamma**2 ) 

#        for j in range(0,len(xin)): # jth peak of this dif
#            
#            if(yin[j] < 0.001):
#                continue
#            
#            delta = np.sqrt( abs( -2.0*sigma**2 *np.log( max(k*sigma*np.sqrt(2.0*np.pi)/(yin[j]*scaling), 0.0000001) ) ) )
#            
#            i_start,i_end = np.searchsorted(xValues,[xin[j]-delta,xin[j]+delta])
#            
#            for i in range(i_start,i_end): # ith x value
#                yout[i] += yin[j] * scaling * temp * np.exp( - (xValues[i] - xin[j])**2 / (2 * sigma**2))

        
        return np.array(yout)
        #TODO: optimize the above method's code to be like the commented-out part at bottom


    def FindCrystalSystem(self):
        """crystal systems from  https://www.geo.arizona.edu/xtal/geos306/geos306-9.htm"""
        #a = self.cell_parameters[0]
        #b = self.cell_parameters[1]
        #c = self.cell_parameters[2]
        #alpha = self.cell_parameters[3]
        #beta = self.cell_parameters[4]
        #gamma = self.cell_parameters[5]
        a,b,c,alpha,beta,gamma = self.cell_parameters_guess
        
        if a == b and a == c and alpha == 90.0 and beta == 90.0 and gamma == 90.0:
            self.crystal_system = "cubic"
            self.num_unique_cell_parameters = 1
        elif a == b and alpha == 90.0 and beta == 90.0 and gamma == 120.0:
            self.crystal_system = "hexagonal"
            self.num_unique_cell_parameters = 2
        elif a == b and a == c and alpha == beta and alpha == gamma:
            self.crystal_system = "rhombohedral"
            self.num_unique_cell_parameters = 2
        elif a == b and alpha == 90.0 and beta == 90.0 and gamma == 90.0:
            self.crystal_system = "tetragonal"
            self.num_unique_cell_parameters = 2
        elif alpha == 90.0 and beta == 90.0 and gamma == 90.0:
            self.crystal_system = "orthorhombic"
            self.num_unique_cell_parameters = 3
        elif beta == 90.0 and gamma == 90.0:
            self.crystal_system = "monoclinic1"
            self.num_unique_cell_parameters = 4
        elif alpha == 90.0 and gamma == 90.0:
            self.crystal_system = "monoclinic2"
            self.num_unique_cell_parameters = 4
        elif beta == 90.0 and gamma == 90.0:
            self.crystal_system = "monoclinic3"
            self.num_unique_cell_parameters = 4
        else:
            self.crystal_system = "triclinic"
            self.num_unique_cell_parameters = 6
    
    def OutputCellParameters(self):
        a = self.cell_parameters_optimized[0]
        b = self.cell_parameters_optimized[1]
        c = self.cell_parameters_optimized[2]
        alpha = self.cell_parameters_optimized[3]
        beta = self.cell_parameters_optimized[4]
        gamma = self.cell_parameters_optimized[5]
        
        if self.crystal_system == "cubic":
            return [a]
        elif self.crystal_system == "hexagonal":
            return [a,c]
        elif self.crystal_system == "rhombohedral":
            return [a,alpha]
        elif self.crystal_system == "tetragonal":
            return [a,c]
        elif self.crystal_system == "orthorhombic":
            return [a,b,c]
        elif self.crystal_system == "monoclinic1":
            return [a,b,c,alpha]
        elif self.crystal_system == "monoclinic2":
            return [a,b,c,beta]
        elif self.crystal_system == "monoclinic3":
            return [a,b,c,gamma]
        else: #self.crystal_system == "triclinic"
            return [a,b,c,alpha,beta,gamma]
        
    def OutputCellParameterBounds(self):
        a = self.cell_parameter_bounds[0]
        b = self.cell_parameter_bounds[1]
        c = self.cell_parameter_bounds[2]
        alpha = self.cell_parameter_bounds[3]
        beta = self.cell_parameter_bounds[4]
        gamma = self.cell_parameter_bounds[5]
        
        if self.crystal_system == "cubic":
            return [a]
        elif self.crystal_system == "hexagonal":
            return [a,c]
        elif self.crystal_system == "rhombohedral":
            return [a,alpha]
        elif self.crystal_system == "tetragonal":
            return [a,c]
        elif self.crystal_system == "orthorhombic":
            return [a,b,c]
        elif self.crystal_system == "monoclinic1":
            return [a,b,c,alpha]
        elif self.crystal_system == "monoclinic2":
            return [a,b,c,beta]
        elif self.crystal_system == "monoclinic3":
            return [a,b,c,gamma]
        else: #self.crystal_system == "triclinic"
            return [a,b,c,alpha,beta,gamma]
        
    def InputCellParameters(self,cp):
        a = None
        b = None
        c = None
        alpha = None
        beta = None
        gamma = None
        
        if self.crystal_system == "cubic":
            a = cp[0]
            b = a
            c = a
            alpha = 90.0
            beta = 90.0
            gamma = 90.0
        elif self.crystal_system == "hexagonal":
            a = cp[0]
            b = a
            c = cp[1]
            alpha = 90.0
            beta = 90.0
            gamma = 120.0
        elif self.crystal_system == "rhombohedral":
            a = cp[0]
            b = a
            c = a
            alpha = cp[1]
            beta = alpha
            gamma = alpha
        elif self.crystal_system == "tetragonal":
            a = cp[0]
            b = a
            c = cp[1]
            alpha = 90.0
            beta = 90.0
            gamma = 90.0            
        elif self.crystal_system == "orthorhombic":
            a = cp[0]
            b = cp[1]
            c = cp[2]
            alpha = 90.0
            beta = 90.0
            gamma = 90.0
        elif self.crystal_system == "monoclinic1":
            a = cp[0]
            b = cp[1]
            c = cp[2]
            alpha = cp[3]
            beta = 90.0
            gamma = 90.0
        elif self.crystal_system == "monoclinic2":
            a = cp[0]
            b = cp[1]
            c = cp[2]
            alpha = 90.0
            beta = cp[3]
            gamma = 90.0
        elif self.crystal_system == "monoclinic3":
            a = cp[0]
            b = cp[1]
            c = cp[2]
            alpha = 90.0
            beta = 90.0
            gamma = cp[3]
        else: #self.crystal_system == "triclinic"
            a = cp[0]
            b = cp[1]
            c = cp[2]
            alpha = cp[3]
            beta = cp[4]
            gamma = cp[5]
        self.cell_parameters_optimized = [a,b,c,alpha,beta,gamma]
    
    def Find2ThetasFromCellParameters(self,true_for_optimized_false_for_guess=True):
        cos = np.cos
        
        if true_for_optimized_false_for_guess:    
            a,b,c,alpha,beta,gamma = self.cell_parameters_optimized
        else:
            a,b,c,alpha,beta,gamma = self.cell_parameters_guess
        
        alpha,beta,gamma = np.deg2rad(alpha), np.deg2rad(beta), np.deg2rad(gamma)
        
        #1. construct G matrix
        G = np.zeros([3,3])
        G = np.asarray([[ a*a,             a*b*cos(gamma),  a*c*cos(beta)  ],
                        [ a*b*cos(gamma),  b*b,             b*c*cos(alpha) ],
                        [ a*c*cos(beta),   b*c*cos(alpha),  c*c            ]])
        G_inv = np.linalg.inv(G)
        
        #2. find d-spacings from G inverse and hkl vectors
        d_spacings = []
        for hkl in self.HKLs:
            d = (hkl.T.dot(G_inv).dot(hkl))**(-0.5)
            d_spacings.append(d)
        
        #3. find 2-thetas from d-spacings
        TwoThetas = [2.0 * np.rad2deg(np.arcsin(self.x_ray_wavelength / (2.0*d))) for d in d_spacings]
        
        self.peak_data[0] = np.asarray(TwoThetas)
    
    def GetUnitCellVolume(self,true_for_optimized_false_for_guess=True):
        if true_for_optimized_false_for_guess:    
            a,b,c,alpha,beta,gamma = self.cell_parameters_optimized
        else:
            a,b,c,alpha,beta,gamma = self.cell_parameters_guess
            
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
        return a*b*c * np.sqrt( 1.0 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) )
    
    
    def ReadFile(self,file_path,print_warnings=True):
        f = open(file_path, "r")
        RIR = 0.0
        data = []
        d_spacings = []
        HKLs = []
        continueReading = True
        readingPeaks = False
        readingChemistry = False
        RIR_found, data_found, chemistry_found = False, False, False
        cell_parameters_found, space_group_found, x_ray_wavelength_found = False, False, False
        d_spacings_found, HKLs_found = False, False
        chemistry = {}
    
        while(continueReading):
            line = f.readline()
            if(line == ""): # end of file
                continueReading = False
            elif(readingPeaks):
                if("2-THETA" in line): # header line, so skip
                    continue
                elif("=" in line):
                    readingPeaks = False
                else:
                    linesplit = line.lstrip().split()
                    x = float(linesplit[0])
                    y = float(linesplit[1])
                    if y != 0.0:
                        data.append([x, y])
                        if len(linesplit) > 3:
                            d_spacings.append(float(linesplit[2]))
                            d_spacings_found = True
                        if len(linesplit) > 5:
                            HKLs.append([int(linesplit[3]),int(linesplit[4]),int(linesplit[5])])
                            HKLs_found = True
                    
            elif(readingChemistry):
                if(line.lstrip() == "" or "#" in line):
                    readingChemistry = False
                else:
                    #for c in chemistryCharactersToDelete:
                    #    line=line.replace(c,"")
                    name,percentage = line.lstrip().split()
                    chemistry[name] = float(percentage)
            elif(line.lstrip().startswith("RIR:")):
                RIR_found = True
                RIR = float(line.lstrip().split()[1])
            elif(line.lstrip().startswith("CELL PARAMETERS:")):
                cell_parameters_found = True
                #dif_cell_parameters.append([float(k) for k in line.lstrip().split()[2:9]])
                self.cell_parameters_guess = [float(k) for k in line.lstrip().split()[2:9]]
                self.cell_parameters_optimized = self.cell_parameters_guess
            elif(line.lstrip().startswith("X-RAY WAVELENGTH:")):
                x_ray_wavelength_found = True
                #dif_x_ray_wavelengths.append(float(line.lstrip().split()[2]))
                self.x_ray_wavelength = float(line.lstrip().split()[2])
            elif("SPACE GROUP:" in line):
                space_group_found = True
                #dif_space_groups.append(line.lstrip().split()[2])
                self.space_group = line.lstrip().split()[-1]
            elif(line.lstrip().startswith("# CHEMISTRY")):
                chemistry_found = True
                readingChemistry = True
            elif(line.lstrip().startswith("# DIF") or line.lstrip().startswith("2-THETA")):
                data_found = True
                readingPeaks = True
        f.close()
        
        if not data_found:
            raise SyntaxError("No data found in dif file: '"+file_path+
                              "'\nFile must contain string '2-THETA' in the line immediately before data.")

        
        data = np.array(data).T
        #data = data.T
        #difXYs.append(data)
        #dif_RIRs.append(RIR)
        #dif_chemistrys.append(chemistry)
        self.peak_data = data
        self.d_spacings = d_spacings
        self.HKLs = np.asarray(HKLs)
        self.RIR = RIR
        self.chemistry = chemistry
        if print_warnings:
            if(not RIR_found):
                print("Warning: no RIR found in dif file: "+file_path)
            if(not data_found):
                print("Warning: no data found in dif file: "+file_path)
                print("Dif data must be preceded by #DIF")
            if(not d_spacings_found):
                print("Warning: no d spacings found in dif file: "+file_path)
            if(not HKLs_found):
                print("Warning: no HKLs found in dif file: "+file_path)
            if(not chemistry_found):
                print("Warning: no chemistry found in dif file: "+file_path)
            if(not cell_parameters_found):
                print("Warning: no cell parameters found in dif file: "+file_path)
            if(not space_group_found):
                print("Warning: no space group found in dif file: "+file_path)
            if(not x_ray_wavelength_found):
                print("Warning: no x-ray wavelength found in dif file: "+file_path)




#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################



class TopLevel(object):
    '''Top level class that wraps everything together. This the only thing the end user interacts with.'''
    
    def __init__(self, input_file_path, twotheta_ranges=[(0.0,90.0)], take_every_nth_point=1, print_warnings = True):
        
        self.Set2ThetaRanges(twotheta_ranges)
        self.twotheta_offset_guess = 0.0
        self.twotheta_offset_optimized = 0.0
        self.refine_twotheta_offset = False
        self.twotheta_offset_bounds = (None,None);
        
        self.print_warnings = print_warnings
        
        self.input_profile = XrayProfile(input_file_path, is_input_profile=True, take_every_nth_point=take_every_nth_point, twotheta_ranges=twotheta_ranges, print_warnings=self.print_warnings)
#        self.input_profile = XrayProfile(input_file_path, is_input_profile=True, take_every_nth_point=take_every_nth_point, twotheta_ranges=twotheta_ranges, twotheta_offset=twotheta_offset, print_warnings=self.print_warnings)
#        self.input_profile.SetScalingGuess(1.0)
        
        self.background_guess = [min(self.input_profile.xy_data[1])]
        self.background_optimized = self.background_guess
        self.background_bounds = [(None,None)]
        self.refine_background = True
        
        self.difs = []
        self.profiles = []
        self.True_for_area_percent_False_for_wt_percent = False
        self.algorithm = "2-step"
        
        self.all_oxide_names = None
        self.calculated_oxide_percentages = None
        
        self.algorithm_run_time = None
        self.algorithm_number_of_iterations = None
        self.norm_between_optimized_profile_and_input = None
        
    def AddPhase(self,file_path,scaling_guess=0.1,FWHM_guess=0.1,FWHM_bounds=(0.05,8.0)):
        #DEPRECATED
        
        added_phase_bool = False
        
        if file_path.endswith(".dif"):
            self.AddDif(file_path,scaling_guess,FWHM_guess,FWHM_bounds)
            added_phase_bool = True
        elif file_path.endswith(".xy") or file_path.endswith(".mdi"):
            self.AddProfile(file_path,scaling_guess)
            added_phase_bool = True
            
        else:
            f_contents = open(file_path, "r").read()
            for line in f_contents.strip().split('\n'):
                if "MAX. ABS. INTENSITY / VOLUME**2:" in line or "# DIF" in line:
                    self.AddDif(file_path,scaling_guess,FWHM_guess,FWHM_bounds)
                    added_phase_bool = True
                    break
                elif "# XY" in line or "# MDI" in line:
                    self.AddProfile(file_path,scaling_guess)
                    added_phase_bool = True
                    break
        
        if not added_phase_bool:
            print("ERROR: no data could be found in file: "+file_path)
            print("""Please indicate the start of XY data or MDI data with #XY or #MDI,
            and for difs please make sure either the file extension is '.dif'
            or it contains the line 'MAX. ABS. INTENSITY / VOLUME**2: ...'""")
    
    def AddDif(self,file_path,use_default=None,scaling=None,scaling_bounds=(0.0,None),refine_scaling=True,FWHM=None,FWHM_bounds=None,refine_FWHM=True,cell_parameters=None,cell_parameter_bounds=None,refine_cell_parameters=True,pseudovoigt_parameter=None,pseudovoigt_parameter_bounds=(0.0,1.0),refine_pseudovoigt_parameter=False):
        d = XrayDif(file_path,print_warnings=self.print_warnings)
        
        if use_default == "crystalline":
            if scaling is None:
                scaling = 0.1
            if FWHM is None:
                FWHM = 0.24
            if FWHM_bounds is None:
                FWHM_bounds = (0.1,0.5)
        elif use_default == "amorphous":
            if scaling is None:
                scaling = 0.3
            if FWHM is None:
                FWHM = 7.0
            if FWHM_bounds is None:
                FWHM_bounds = (2.0,20.0)
        else:
            # assume crystalline but have the peak width bounds allow it to become amorphous
            if scaling is None:
                scaling = 0.1
            if FWHM is None:
                FWHM = 0.24
            if FWHM_bounds is None:
                FWHM_bounds = (0.1,20.0)
        
        if pseudovoigt_parameter == None:
            pseudovoigt_parameter = 0.0
        
        if FWHM_bounds[0] < 0.01:
            FWHM_bounds = (0.01,FWHM_bounds[1])
            
        d.SetScalingGuess(scaling)
        d.SetScalingBounds(scaling_bounds)
        d.refine_scaling = refine_scaling
        d.SetFWHMGuess(FWHM)
        d.SetFWHMBounds(FWHM_bounds)
        d.refine_FWHM = refine_FWHM
        d.SetPseudoVoigtParameterGuess(pseudovoigt_parameter)
        d.SetPseudoVoigtParameterBounds(pseudovoigt_parameter_bounds)
        d.refine_pseudovoigt_parameter = refine_pseudovoigt_parameter
        
        d.cell_parameters_original = list(d.cell_parameters_guess) #for copypasta
        if cell_parameters is not None:
            d.SetCellParametersGuess(cell_parameters)
            d.SetCellParametersOptimized(cell_parameters)
            d.Find2ThetasFromCellParameters()
        if cell_parameter_bounds is not None:
            d.SetCellParameterBoundsUsingValues(cell_parameter_bounds)
        else:
            d.SetCellParameterBoundsUsingPercentage(10.0)
        d.refine_cell_parameters = refine_cell_parameters
        
        self.difs.append(d)
    
    def AddProfile(self,file_path,scaling=0.1,scaling_bounds=(0.0,None),refine_scaling=True,twotheta_offset=0.0):
        p = XrayProfile(file_path,is_input_profile=False,twotheta_offset=0.0,print_warnings=self.print_warnings)
        p.SetScalingGuess(scaling)
        p.SetScalingBounds(scaling_bounds)
        p.refine_scaling = refine_scaling
        p.twotheta_offset = twotheta_offset
        self.profiles.append(p)
    
    def SetBackgroundGuess(self, background):
        self.background_guess = background
        self.background_optimized = background
        self.background_bounds = [(None,None)] * len(background)
    
    def SetBackgroundBounds(self,background_bounds):
        self.background_bounds = background_bounds
    
    def RefineBackground(self, boolean):
        self.refine_background = boolean
        
    def Set2ThetaOffsetGuess(self, offset):
        self.twotheta_offset_guess = offset
        self.twotheta_offset_optimized = offset
    
    def Set2ThetaOffsetBounds(self,offset_bounds):
        self.twotheta_offset_bounds = offset_bounds

    def Refine2ThetaOffset(self,boolean):
        self.refine_twotheta_offset = boolean
    
    def Set2ThetaRanges(self, list_of_ranges):
        # Can input either like [start,end] or [[start1,end1],(start2,end2), ...]
        # but will always be stored like [(start1,end1),(start2,end2), ...]
        if list_of_ranges is None:
            self.twotheta_ranges = [(0.0,90.0)]
            return
        if isinstance(list_of_ranges[0], (int, long, float)):
            # list_of_ranges was input like [start,end]
            self.twotheta_ranges = [(float(list_of_ranges[0]),float(list_of_ranges[1]))]
            return
        else:
            # list_of_ranges was input like [[start1,end1],...]
            self.twotheta_ranges = []
            for pair in list_of_ranges:
                self.twotheta_ranges += [(float(pair[0]),float(pair[1]))]
            return
            
    
    def SumDifsAndProfiles(self,true_for_optimized_false_for_guess=True,true_for_chopped_data_false_for_full_original=True):
        
        if true_for_chopped_data_false_for_full_original:
            inputXY = self.input_profile.xy_data
        else:
            inputXY = self.input_profile.xy_data_for_plotting
        
        if true_for_optimized_false_for_guess:
            background = self.background_optimized
        else:
            background = self.background_guess
        
#        y_out = np.full_like(inputXY[0],background[0])
#        for degree in range(len(background))[1:]:
#            y_out = [y_out[i] + background[degree]*self.input_profile.xy_data[0]**degree for i in range(len(y_out))]
        y_out = np.array(map(lambda x: sum([background[degree]*x**degree for degree in range(len(background))]), list(inputXY[0])))
        
        for d in self.difs:
            y_out = y_out + d.GetBroadenedPattern(inputXY[0],true_for_optimized_false_for_guess=true_for_optimized_false_for_guess)
        for p in self.profiles:
            y_out = y_out + p.GetScaledPattern(inputXY[0],true_for_optimized_false_for_guess=true_for_optimized_false_for_guess)

        return y_out

    def InitializeBackgroundToMinimumOfInputPattern(self, degree = 0):
        self.background_guess = [min(self.input_profile.xy_data[1])]
        for i in range(degree):
            self.background_guess += [0.0]
        self.background_optimized = self.background_guess

    def GetNorm(self):
        y_out = self.SumDifsAndProfiles()
        self.norm_between_optimized_profile_and_input = np.linalg.norm(self.input_profile.xy_data[1] - y_out)
        return self.norm_between_optimized_profile_and_input
    
    
#    def CalculateChemistry(self):
#        #find all_oxides
#        
#        profiles = self.profiles
#        difs = self.difs
#        input_file_chemistry = self.input_profile.chemistry
#        all_oxide_names = list(set( key for dic in [p.chemistry for p in profiles]+[d.chemistry for d in difs]+[input_file_chemistry] for key in dic.keys()))
#        self.all_oxide_names = all_oxide_names
#        start_2theta = min(self.input_profile.xy_data[0])
#        end_2theta = max(self.input_profile.xy_data[0])
#        
#        oxide_sums = {}
#        for oxide in all_oxide_names:
#            calculated_absolute_amount = 0.0
#            for p in profiles:
#                if(p.RIR > 0.0):
#                    calculated_absolute_amount += p.GetArea(self.input_profile.xy_data[0]) * p.chemistry.get(oxide,0.0) / p.RIR
#            for d in difs:
#                if(d.RIR > 0.0):
#                    calculated_absolute_amount += d.GetArea(start_2theta = start_2theta,end_2theta = end_2theta) * d.chemistry.get(oxide,0.0) / d.RIR
#            oxide_sums[oxide] = calculated_absolute_amount
#        sum_of_oxide_sums = sum(oxide_sums.values())
#        
#        for oxide in all_oxide_names:
#            target_percentage = input_file_chemistry.get(oxide,0.0)
#            calculated_percentage = oxide_sums[oxide] / sum_of_oxide_sums * 100.0
        
        
    def CalculateChemistry(self):
        profiles = self.profiles
        difs = self.difs
        input_file_chemistry = self.input_profile.chemistry
        all_oxide_names = list(set( key for dic in [p.chemistry for p in profiles]+[d.chemistry for d in difs]+[input_file_chemistry] for key in dic.keys()))
        self.all_oxide_names = all_oxide_names
        start_2theta = min(self.input_profile.xy_data[0])
        end_2theta = max(self.input_profile.xy_data[0])
        
        calculated_oxide_percentages = {}
        
        #calculate sum over all phases of area/RIR
        sum_of_areas = 0.0
        for p in [p for p in profiles if p.RIR != 0.0]:
#            p.SetAbsoluteContribution( p.GetArea(self.input_profile.xy_data[0])/p.RIR )
            p.SetAbsoluteContribution( p.scaling_optimized*100.0/p.RIR )
            sum_of_areas += p.GetAbsoluteContribution()
        for d in [d for d in difs if d.RIR != 0.0]:
#            d.SetAbsoluteContribution( d.GetArea(start_2theta = start_2theta, end_2theta = end_2theta)/d.RIR )
            d.SetAbsoluteContribution( d.scaling_optimized*100.0/d.RIR )
            sum_of_areas += d.GetAbsoluteContribution()
        
        phases = [p for p in profiles if p.RIR != 0.0] + [d for d in difs if d.RIR != 0.0]
        
        #calculate relative contributions
        for p in phases:
            p.SetRelativeContribution( 100.0 * p.GetAbsoluteContribution() / sum_of_areas )
        
        for oxide in all_oxide_names:
            calculated_percentage = 0.0
            for p in phases:
                calculated_percentage += p.GetRelativeContribution() * p.chemistry.get(oxide,0.0) / 100.0
            calculated_oxide_percentages[oxide] = calculated_percentage
            
        self.calculated_oxide_percentages = calculated_oxide_percentages

    def OutputCSV(self,file_path="output - all phases as profiles.csv",only_optimized_xy=False):

        inputXY = self.input_profile.xy_data_for_plotting
        
        optimizedY = self.SumDifsAndProfiles(true_for_optimized_false_for_guess=True,true_for_chopped_data_false_for_full_original=False)
        

        if only_optimized_xy:
            f = open(file_path,"w+")
            x_values = inputXY[0]
            y_values = optimizedY
            for i in range(len(x_values)):
                f.write(str(x_values[i])+","+str(y_values[i])+"\n")
            f.close()
        
        else:
            background = np.array(map(lambda x: sum([self.background_optimized[degree]*x**degree for degree in range(len(self.background_optimized))]), list(inputXY[0])))
            x_values = inputXY[0]
            y_values = [inputXY[1],optimizedY,background]
            header_line = 'input_x,input_y,optimized_y,background'
            for d in self.difs:
                y_values += [d.GetBroadenedPattern(inputXY[0])]
                header_line += ',' + d.file_nickname
            for p in self.profiles:
                y_values += [p.GetScaledPattern(inputXY[0])]
                header_line += ',' + p.file_nickname
            
            f = open(file_path,"w+")
            f.write(header_line+'\n')
            for i in range(len(x_values)):
                f.write(str(x_values[i]))
                for j in range(len(y_values)):
                    f.write(','+str(y_values[j][i]))
                f.write('\n')
            f.close()
            
    def DoFinalPlot(self,range_to_plot=None,linewidth=0.75):
        
        self.CalculateChemistry()
        
        difs = self.difs
        profiles = self.profiles
        
        my_plot = plt.figure(1,figsize=(30, 15))
        
        inputXY = np.array(self.input_profile.xy_data_for_plotting)
        if self.twotheta_offset_optimized != 0.0:
            inputXY[0] = [inputXY[0][i] + self.twotheta_offset_optimized for i in range(len(inputXY[0]))]
            for p in self.profiles:
                p.GetInterpolatedXYData(inputXY[0])
        
        if range_to_plot is None:
            if self.twotheta_ranges == [(0.0,90.0)]:
                x_min, x_max = min(inputXY[0]), max(inputXY[0])
                global_multiplier = 100.0/max(inputXY[1])
            else:
                x_min = min([pair[0] for pair in self.twotheta_ranges])
                x_max = max([pair[1] for pair in self.twotheta_ranges])
                global_multiplier = 100.0/max([inputXY[1][i] for i in range(len(inputXY[0])) if (inputXY[0][i] >= x_min and inputXY[0][i] <= x_max)])
        else:
            x_min, x_max = range_to_plot[0], range_to_plot[1]
            global_multiplier = 100.0/max([inputXY[1][i] for i in range(len(inputXY[0])) if (inputXY[0][i] >= x_min and inputXY[0][i] <= x_max)])
        
        optimizedY = self.SumDifsAndProfiles(true_for_optimized_false_for_guess=True,true_for_chopped_data_false_for_full_original=False)
        
#        background = np.full_like(inputXY[1],self.background_optimized[0])
        background = np.array(map(lambda x: sum([self.background_optimized[degree]*x**degree for degree in range(len(self.background_optimized))]), list(inputXY[0])))

        plt.subplot(222)
        for p in profiles:
            if(p.RIR < 0.001):
                plt.plot(inputXY[0],(p.GetScaledPattern(inputXY[0])+background)*global_multiplier,label=p.file_nickname,linewidth=linewidth)
            else:
                if(self.True_for_area_percent_False_for_wt_percent):
                    plt.plot(inputXY[0],(p.GetScaledPattern(inputXY[0])+background)*global_multiplier,label=p.file_nickname+" "+
                             "{0:.1f}".format(p.GetRelativeContribution())+"%",linewidth=linewidth)
                else:
#                    plt.plot(inputXY[0],(p.GetScaledPattern(inputXY[0])+background)*global_multiplier,label=p.file_nickname+" "+
#                             "{0:.1f}".format(p.GetRelativeContribution())+" wt%",linewidth=linewidth)
                    plt.plot(inputXY[0],(p.GetScaledPattern(inputXY[0])+background)*global_multiplier,label=p.file_nickname,linewidth=linewidth)
    
        for d in difs:
            if(d.RIR < 0.001):
                plt.plot(inputXY[0],(d.GetBroadenedPattern(inputXY[0])+background)*global_multiplier,label=d.file_nickname,linewidth=linewidth)
            else:
                if(self.True_for_area_percent_False_for_wt_percent):
                    plt.plot(inputXY[0],(d.GetBroadenedPattern(inputXY[0])+background)*global_multiplier,label=d.file_nickname+" "+
                        "{0:.1f}".format(d.GetRelativeContribution())+"%",linewidth=linewidth)
                else:
#                    plt.plot(inputXY[0],d.GetBroadenedPattern(inputXY[0])+background,label=d.file_nickname+" "+
#                        "{0:.2f}".format(d.GetRelativeContribution())+" wt%",linewidth=linewidth)
                    plt.plot(inputXY[0],(d.GetBroadenedPattern(inputXY[0])+background)*global_multiplier,label=d.file_nickname,linewidth=linewidth)
        plt.plot(inputXY[0],background*global_multiplier,label="background",linestyle="dotted",color="red",linewidth=linewidth)
        plt.plot(inputXY[0],inputXY[1]*global_multiplier,label="input:"+self.input_profile.file_nickname,color="black",linewidth=linewidth)
        plt.plot(inputXY[0],optimizedY*global_multiplier,label="optimized profile",color="magenta",linewidth=linewidth)
        plt.grid(True)
        plt.legend(bbox_to_anchor=(1.0, -0.05))
        plt.xlabel('2-theta (deg)')
        plt.ylabel('intensity')
        difference = (inputXY[1] - optimizedY)*global_multiplier
        offsetline = np.full_like(inputXY[0],110.0)
        plt.plot(inputXY[0], difference+offsetline, linestyle="solid", color="pink",linewidth=linewidth)
        plt.plot(inputXY[0], offsetline, linestyle="dotted", color="red")
        plt.axis([x_min, x_max, 0.0 ,120.0])
        plt.show()

    def DoGuessPlot(self,range_to_plot=None,linewidth=0.75):
        
        difs = self.difs
        profiles = self.profiles
        
        inputXY = np.array(self.input_profile.xy_data_for_plotting)
        if self.twotheta_offset_optimized != 0.0:
            inputXY[0] = [inputXY[0][i] + self.twotheta_offset_guess for i in range(len(inputXY[0]))]
            for p in self.profiles:
                p.GetInterpolatedXYData(inputXY[0])

        if range_to_plot is None:
            if self.twotheta_ranges == [(0.0,90.0)]:
                x_min, x_max = min(inputXY[0]), max(inputXY[0])
                global_multiplier = 100.0/max(inputXY[1])
            else:
                x_min = min([pair[0] for pair in self.twotheta_ranges])
                x_max = max([pair[1] for pair in self.twotheta_ranges])
                global_multiplier = 100.0/max([inputXY[1][i] for i in range(len(inputXY[0])) if (inputXY[0][i] >= x_min and inputXY[0][i] <= x_max)])
        else:
            x_min, x_max = range_to_plot[0], range_to_plot[1]
            global_multiplier = 100.0/max([inputXY[1][i] for i in range(len(inputXY[0])) if (inputXY[0][i] >= x_min and inputXY[0][i] <= x_max)])

#        background = np.full_like(inputXY[1],self.background_guess[0])
        background = np.array(map(lambda x: sum([self.background_optimized[degree]*x**degree for degree in range(len(self.background_optimized))]), list(inputXY[0])))
        
#        for p in profiles:
#            p.GetInterpolatedXYData(inputXY[0])
            
        guessY = self.SumDifsAndProfiles(true_for_optimized_false_for_guess=False,true_for_chopped_data_false_for_full_original=False)
        
        my_plot = plt.figure(1,figsize=(30, 15))
        plt.subplot(222)
        for p in profiles:
            plt.plot(inputXY[0],(p.GetScaledPattern(inputXY[0],true_for_optimized_false_for_guess=False)+background)*global_multiplier,label=p.file_nickname,linewidth=linewidth)
        for d in difs:
            plt.plot(inputXY[0],(d.GetBroadenedPattern(inputXY[0],true_for_optimized_false_for_guess=False)+background)*global_multiplier,label=d.file_nickname,linewidth=linewidth)
        plt.plot(inputXY[0],background*global_multiplier,label="background",linestyle="dotted",color="red",linewidth=linewidth)
        plt.plot(inputXY[0],inputXY[1]*global_multiplier,label="input:"+self.input_profile.file_nickname,color="black",linewidth=linewidth)
        plt.plot(inputXY[0],guessY*global_multiplier,label="guess profile",color="magenta",linewidth=linewidth)
#        if range_to_plot is None:
#            if self.twotheta_ranges == [(0.0,90.0)]:
#                plt.axis([min(inputXY[0]), max(inputXY[0]), 0.0 ,120.0])
#            else:
#                x_min = min([pair[0] for pair in self.twotheta_ranges])
#                x_max = max([pair[1] for pair in self.twotheta_ranges])
#                plt.axis([x_min, x_max, 0.0 ,120.0])
#        else:
#            plt.axis([range_to_plot[0], range_to_plot[1], 0.0 ,120.0])
        plt.grid(True)
        plt.legend(bbox_to_anchor=(1.0, -0.05))
        plt.xlabel('2-theta (deg)')
        plt.ylabel('intensity')
        difference = (inputXY[1] - guessY)*global_multiplier
        offsetline = np.full_like(inputXY[0],110.0)
        plt.plot(inputXY[0], difference+offsetline, linestyle="solid", color="pink",linewidth=linewidth)
        plt.plot(inputXY[0], offsetline, linestyle="dotted", color="red")
        plt.axis([x_min, x_max, 0.0 ,120.0])
        plt.show()
    
    def PrintParameterResults(self):
        self.CalculateChemistry()
        print("\nParameter results:")
#        for d in self.difs:
#            print(d.file_path)
#            print('  optimized FWHM         '+str(d.FWHM_optimized))
##            if d.pseudovoigt_parameter_optimized != 0.0:
#            print('  opt pseudovoigt param  '+str(d.pseudovoigt_parameter_optimized))
#            print('  optimized scaling      '+str(d.scaling_optimized))
#            print('  optimized cell params  '+str(d.cell_parameters_optimized))
#            if d.relative_contribution is None:
#                print('  relative contribution  None')
#            else:
#                print('  relative contribution  '+str(d.relative_contribution)+' wt%')
#        for p in self.profiles:
#            print(p.file_path)
#            print('  optimized scaling      '+str(p.scaling_optimized))
#            if p.relative_contribution is None:
#                print('  relative contribution  None')
#            else:
#                print('  relative contribution  '+str(p.relative_contribution)+' wt%')
#        print('optimized background     '+str(self.background_optimized[0]))
        for d in self.difs:
            print(d.file_path)
            print('  FWHM                   '+'{0:.4f}'.format(d.FWHM_optimized))
            print('  pseudovoigt parameter  '+'{0:.4f}'.format(d.pseudovoigt_parameter_optimized).replace('0.0000','0.0'))
            print('  scaling                '+'{0:.4f}'.format(d.scaling_optimized))
#            print('  opt cell parameters    '+str(["{0:0.4f}".format(i) for i in d.cell_parameters_optimized]).replace("'", "").replace('0.0000','0').replace(',',' '))
            print('  cell parameters        '+str(["{0:0.4f}".format(i) for i in d.cell_parameters_optimized]).replace("'", "").replace('0.0000','0').replace(', ',', ,').replace(']',', ]'))
            print('  unit cell volume       '+'{0:.4f}'.format(d.GetUnitCellVolume()))
            if d.relative_contribution is None:
                print('  relative contribution  None')
            else:
                print('  relative contribution  '+'{0:.4f}'.format(d.relative_contribution)+' wt%')
        for p in self.profiles:
            print(p.file_path)
            print('  scaling                '+'{0:.4f}'.format(p.scaling_optimized))
            if p.relative_contribution is None:
                print('  relative contribution  None')
            else:
                print('  relative contribution  '+'{0:.4f}'.format(p.relative_contribution)+' wt%')
                
#        print('optimized background     '+'{0:.4f}'.format(self.background_optimized[0]))
        print_string = 'background               '+'{0:.4f}'.format(self.background_optimized[0])
        if self.refine_twotheta_offset:
            print('2theta offset            '+'{0:.4f}'.format(self.twotheta_offset_optimized))
        for degree in range(len(self.background_optimized))[1:]:
            print_string += ' + '+'{0:.8f}'.format(self.background_optimized[degree])+' x^'+str(degree)
        print(print_string)
        
    def PrintParameterResultsTable(self):
        self.CalculateChemistry()
        print('\nname,wt%,a,err a,b,err b,c,err c,alpha,err alpha,beta,err beta,gamma,err gamma,vol,err vol,FWHM')
        for d in self.difs:
            print(d.file_nickname+','+str(d.relative_contribution)+','+str(d.cell_parameters_optimized[0])+',,'+str(d.cell_parameters_optimized[1])+',,'+
                  str(d.cell_parameters_optimized[2])+',,'+str(d.cell_parameters_optimized[3])+',,'+str(d.cell_parameters_optimized[4])+',,'+
                  str(d.cell_parameters_optimized[5])+',,'+str(d.GetUnitCellVolume())+',,'+str(d.FWHM_optimized))
        for p in self.profiles:
            print(p.file_nickname+','+str(p.relative_contribution)+',,,,,,,,,,,,,,,')
        
    def PrintChemistryResults(self):
        self.CalculateChemistry()
        print("\nChemistry results:")
#        print("  oxide:".ljust(10)+"calculated chemistry:".ljust(25)+"input file chemistry:".ljust(25)+"difference:".ljust(25))
        print("  "+"oxide:".ljust(12)+"calculated:".ljust(15)+"input file:".ljust(15)+"difference:".ljust(15))
        for oxide in self.all_oxide_names:
            calculated_percent = self.calculated_oxide_percentages.get(oxide,0.0)
            input_percent = self.input_profile.chemistry.get(oxide,0.0)
#            print("  "+ oxide.ljust(10) + (str(calculated_percent)+" wt%").ljust(25) + (str(input_percent)+" wt%").ljust(25)+("%+f" %(calculated_percent-input_percent)+" wt%").ljust(25) )
            print("  "+ oxide.ljust(12) + ('{0:.4f}'.format(calculated_percent)+" wt%").ljust(15) + ('{0:.4f}'.format(input_percent)+" wt%").ljust(15)+('{0:+.4f}'.format(calculated_percent-input_percent)+" wt%").ljust(15) )
    
    def PrintAlgorithmStats(self):
        print("\nAlgorithm stats:")
        print("  algorithm run time    "+str(self.algorithm_run_time))
        print("  number of iterations  "+str(self.algorithm_number_of_iterations))
        print("  norm with input data  "+str(self.GetNorm()))
    
    def PrintCopyPasta(self):
        print("\nnice copypasta:")
        for d in self.difs:
            print_string = "apc.AddDif('"+d.file_path+"'"
            print_string += ",scaling="+str(d.scaling_optimized)
            if d.scaling_bounds != (0.0,None):
                print_string += ",scaling_bounds="+str(d.scaling_bounds)
            if d.refine_scaling == False:
                print_string += ",refine_scaling=False"
            print_string += ",FWHM="+str(d.FWHM_optimized)
            print_string += ",FWHM_bounds="+str(d.FWHM_bounds)
            if d.refine_FWHM == False:
                print_string += ",refine_FWHM=False"
            if d.cell_parameters_optimized != d.cell_parameters_original:
                print_string += ",cell_parameters="+str(d.cell_parameters_optimized)
#            if d.cell_parameter_bounds != [(0.1,None),(0.1,None),(0.1,None),(30.0,None),(30.0,None),(30.0,None)]:
#            if d.cell_parameter_bounds != [(d.cell_parameters_optimized[i]*(1.0-percentage/100.0),self.cell_parameters_optimized[i]*(1.0+percentage/100.0)) for i in range(6)]
            if d.cell_parameter_bounds != [(d.cell_parameters_optimized[i]*0.90,d.cell_parameters_optimized[i]*1.1) for i in range(6)]:
                print_string += ",cell_parameter_bounds="+str(d.cell_parameter_bounds)
            if d.refine_cell_parameters == False:
                print_string += ",refine_cell_parameters=False"
            if d.pseudovoigt_parameter_optimized != 0.0:
                print_string += ",pseudovoigt_parameter="+str(d.pseudovoigt_parameter_optimized)
            if d.pseudovoigt_parameter_bounds != (0.0,1.0):
                print_string += ",pseudovoigt_parameter_bounds="+str(d.pseudovoigt_parameter_bounds)
            if d.refine_pseudovoigt_parameter == True:
                print_string += ",refine_pseudovoigt_parameter=True"
            print_string += ")"
            print(print_string)
#            print("apc.AddDif('"+d.file_path+"',scaling_guess="+str(d.scaling_optimized)+",FWHM_guess="+str(d.FWHM_optimized)+",FWHM_bounds="+str(d.FWHM_bounds)+",cell_parameters_guess="+str(d.cell_parameters_optimized)+",refine_cell_parameters="+str(d.refine_cell_parameters)+")")
        for p in self.profiles:
            print_string = "apc.AddProfile('"+p.file_path+"'"
            print_string += ",scaling="+str(p.scaling_optimized)
            if p.scaling_bounds != (0.0,None):
                print_string += ",scaling_bounds="+str(p.scaling_bounds)
            if p.refine_scaling == False:
                print_string += ",refine_scaling=False"
            if p.twotheta_offset != 0.0:
                print_string += ",twotheta_offset="+str(p.twotheta_offset)
            print_string += ")"
            print(print_string)
#            print("apc.AddProfile('"+p.file_path+"',scaling_guess="+str(p.scaling_optimized)+")")
        print("apc.SetBackgroundGuess("+str(self.background_optimized)+")")
        if not self.refine_background:
            print("apc.RefineBackground(False)")
        if self.twotheta_offset_optimized != 0.0:
            print("apc.Set2ThetaOffset("+str(self.twotheta_offset_optimized)+")")
        if self.refine_twotheta_offset:
            print("apc.Refine2ThetaOffset(True)")
    
    def PrintEverything(self):
        self.PrintAlgorithmStats()
        self.PrintParameterResults()
        self.PrintParameterResultsTable()
        self.PrintChemistryResults()
        self.PrintCopyPasta()
        self.OutputCSV()
    
    
    def DoOptimization(self,algorithm="2-step",RIR_implementation="full integrated area",chemistry_constraint_a=None,chemistry_constraint_b=None,print_start_and_finish=True):
        """possible_algorithms:
            "2-step",
            "all nonlinear",
        """
        possible_RIR_implementations = [
                "full integrated area",
                "area of highest peak only"
                ]
        possible_algorithms = [
                "2-step",
                "all nonlinear",
                ]
        
        if algorithm not in possible_algorithms:
            print("ERROR: not a valid algorithm name. Please use one of the following:")
            print(possible_algorithms)
            return
        else:
            self.algorithm = algorithm
        
        if RIR_implementation not in possible_RIR_implementations:
            print("ERROR: not a valid RIR implementation. Please use one of the following:")
            print(possible_RIR_implementations)
            return
        else:
            self.RIR_implementation = RIR_implementation

        if print_start_and_finish:
            print("started optimization")
        
        self.RunFlexibleOptimization(algorithm,chemistry_constraint_a,chemistry_constraint_b)
        
        if print_start_and_finish:
            print("finished optimization")
    
    def RunFlexibleOptimization(self, algorithm="2-step", chemistry_constraint_a=None, chemistry_constraint_b=None):
        
        # refine_scaling not implemented for 2-step yet
        # profiles only not implemented for 2-step yet
        
        inputXY = np.array(self.input_profile.xy_data)
        
        x = []
        bounds = []
        
        if self.refine_twotheta_offset:
            x += [self.twotheta_offset_optimized]
            bounds += [self.twotheta_offset_bounds]
        if self.twotheta_offset_optimized != 0.0:
            inputXY[0] = [inputXY[0][i] + self.twotheta_offset_optimized for i in range(len(inputXY[0]))]
        
        if algorithm == "all nonlinear":
            if self.refine_background:
                x += self.background_optimized
                bounds += self.background_bounds

        #FIXME: do interpolation earlier
        for p in self.profiles:
            p.GetInterpolatedXYData(inputXY[0])
            if algorithm == "all nonlinear":
                if p.refine_scaling:
                    x += [p.scaling_optimized]
                    bounds += [p.scaling_bounds]
        
        # determine crystal system of each dif
        for d in self.difs:
            d.FindCrystalSystem()
            if algorithm == "all nonlinear":
                if d.refine_scaling:
                    x += [d.scaling_optimized]
                    bounds += [d.scaling_bounds]
            if d.refine_FWHM:
                x += [d.FWHM_optimized]
                bounds += [d.FWHM_bounds]
            if d.refine_pseudovoigt_parameter:
                x += [d.pseudovoigt_parameter_optimized]
                bounds += [d.pseudovoigt_parameter_bounds]
            if d.refine_cell_parameters:
                x += d.OutputCellParameters()
                bounds += d.OutputCellParameterBounds()
        
        # x is [twotheta_offset, background, profile scalings, 
        #       repeated units of {dif scaling, dif FWHM, dif pseudovoigt parameter, dif cell parameters} ]
        #################################################
        # run optimization
        t_start = time()
        if algorithm == "2-step":
            if len(x) > 0:
                sol = minimize( self.TwoStep, x, bounds = bounds, method='L-BFGS-B' )
                self.algorithm_number_of_iterations = sol.nit
            else:
                self.algorithm_number_of_iterations = 0
        if algorithm == "all nonlinear":
            sol = minimize( self.AllNonlinear, x, bounds = bounds, method='L-BFGS-B' )
            self.algorithm_number_of_iterations = sol.nit
        self.algorithm_run_time = time() - t_start
#        self.algorithm_number_of_iterations = sol.nit
        #################################################
        
        # unpack results
        # convert x array to list to be able to pop elements off
        if len(x) > 0:
            xlist = list(sol.x)
        
            if self.refine_twotheta_offset:
                self.twotheta_offset_optimized = xlist.pop(0)
                inputXY = np.array(self.input_profile.xy_data)
                if self.twotheta_offset_optimized != 0.0:
                    inputXY[0] = [inputXY[0][i] + self.twotheta_offset_optimized for i in range(len(inputXY[0]))]
                for p in self.profiles:
                    p.GetInterpolatedXYData(inputXY[0])
            
            if algorithm == "all nonlinear":
                if self.refine_background:
                    self.background_optimized = [xlist.pop(0) for i in range(len(self.background_optimized))]
            
            if algorithm == "all nonlinear":
                for p in self.profiles:
                    if p.refine_scaling:
                        p.scaling_optimized = xlist.pop(0)
            
            for d in self.difs:
                if algorithm == "all nonlinear":
                    if d.refine_scaling:
                        d.scaling_optimized = xlist.pop(0)
                if d.refine_FWHM:
                    d.FWHM_optimized = xlist.pop(0)
                if d.refine_pseudovoigt_parameter:
                    d.pseudovoigt_parameter_optimized = xlist.pop(0)
                if d.refine_cell_parameters:
                    d.InputCellParameters(xlist[0:d.num_unique_cell_parameters])
                    del xlist[0:d.num_unique_cell_parameters]
                d.Find2ThetasFromCellParameters()
        
        if algorithm == "2-step":
            # calculate final scalings
            backgroundYs = [np.array(map(lambda x: x**degree, list(inputXY[0]))) for degree in range(len(self.background_optimized))]
            profileYs = [p.y_data_interpolated for p in self.profiles]
            difYs = [self.difs[i].GetBroadenedPattern(inputXY[0],scaling=1.0) for i in range(len(self.difs))]
            
            A = np.array( backgroundYs + profileYs + difYs ).T
            
            solution = np.linalg.lstsq(A, inputXY[1])
            
            xlstsq = list(solution[0])
        
            self.background_optimized = [xlstsq.pop(0) for i in range(len(self.background_optimized))]
            
            for p in self.profiles:
                p.scaling_optimized = xlstsq.pop(0)
            for d in self.difs:
                d.scaling_optimized = xlstsq.pop(0)
        
    def TwoStep(self, x):
        inputXY = np.array(self.input_profile.xy_data)
        
        xlist = list(x)
        
        if self.refine_twotheta_offset:
            self.twotheta_offset_optimized = xlist.pop(0)
            if self.twotheta_offset_optimized != 0.0:
                inputXY[0] = [inputXY[0][i] + self.twotheta_offset_optimized for i in range(len(inputXY[0]))]
            #TODO: I unindented these two lines below. Should they be indented?
#                for p in self.profiles:
#                    p.GetInterpolatedXYData(inputXY[0])
            for p in self.profiles:
                p.GetInterpolatedXYData(inputXY[0])
        
        for d in self.difs:
            if d.refine_FWHM:
                d.FWHM_optimized = xlist.pop(0)
            if d.refine_pseudovoigt_parameter:
                d.pseudovoigt_parameter_optimized = xlist.pop(0)
            if d.refine_cell_parameters:
                d.InputCellParameters(xlist[0:d.num_unique_cell_parameters])
                del xlist[0:d.num_unique_cell_parameters]
            d.Find2ThetasFromCellParameters()
        
        # calculate final scalings
        backgroundYs = [np.array(map(lambda x: x**degree, list(inputXY[0]))) for degree in range(len(self.background_optimized))]
        profileYs = [p.y_data_interpolated for p in self.profiles]
        difYs = [self.difs[i].GetBroadenedPattern(inputXY[0],scaling=1.0) for i in range(len(self.difs))]
        
        A = np.array( backgroundYs + profileYs + difYs ).T
        
        solution = np.linalg.lstsq(A, inputXY[1])
        
        xsol = list(solution[0])
        
        self.background_optimized = [xsol.pop(0) for i in range(len(self.background_optimized))]
        
        for p in self.profiles:
            p.scaling_optimized = xsol.pop(0)
        for d in self.difs:
            d.scaling_optimized = xsol.pop(0)
        
        y_out = self.SumDifsAndProfiles()
        return np.linalg.norm(inputXY[1] - y_out)

    def AllNonlinear(self, x):
        inputXY = np.array(self.input_profile.xy_data)
        
        xlist = list(x)
        # x is [twotheta_offset, background, profile scalings, 
        #       repeated units of {dif scaling, dif FWHM, dif pseudovoigt parameter, dif cell parameters} ]

        if self.refine_twotheta_offset:
            self.twotheta_offset_optimized = xlist.pop(0)
            if self.twotheta_offset_optimized != 0.0:
                inputXY[0] = [inputXY[0][i] + self.twotheta_offset_optimized for i in range(len(inputXY[0]))]
            for p in self.profiles:
                p.GetInterpolatedXYData(inputXY[0])
        
        if self.refine_background:
            self.background_optimized = [xlist.pop(0) for i in range(len(self.background_optimized))]
    
        for p in self.profiles:
            if p.refine_scaling:
                p.scaling_optimized = xlist.pop(0)
        
        for d in self.difs:
            if d.refine_scaling:
                d.scaling_optimized = xlist.pop(0)
            if d.refine_FWHM:
                d.FWHM_optimized = xlist.pop(0)
            if d.refine_pseudovoigt_parameter:
                d.pseudovoigt_parameter_optimized = xlist.pop(0)
            if d.refine_cell_parameters:
                d.InputCellParameters(xlist[0:d.num_unique_cell_parameters])
                del xlist[0:d.num_unique_cell_parameters]
            d.Find2ThetasFromCellParameters()
        
        y_out = self.SumDifsAndProfiles()
        return np.linalg.norm(inputXY[1] - y_out)
    

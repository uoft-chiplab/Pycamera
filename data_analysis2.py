"""
Important and useful data analysis functions.   Bunch them together in this
file for easy import into some namespace later on.  Plan is to call these
function from pycamera.py in the analysis thread.
"""

from pylab import *
from scipy import *
from matplotlib.mlab import find

import pdb

# --- Test function definition.  -------- 
def marcius():
    print "the Q"
    return "the M"

# --- 1D gaussian fit to 2D image data array ---------
def gfit1D(Xdata, Gdata, p0):
    """ Fit data 'Gdata' to a 1D gaussian over the data range 'Xdata'.
        
        p0 is the array (3 values) of inital guesses
        p0[0] = amplitude
        p0[1] = width**2
        p0[2] = centre 
        p0[3] = constant background
    """
    # 1D gaussian definition
    fitfunc = lambda p, x: p[0]*exp( -(x-p[2])**2 / (2*p[1]) ) + p[3]
    # error function
    errfunc = lambda p, x, val: fitfunc(p,x) - val  
    # perform the fit
    fit,success = optimize.leastsq(errfunc, p0.copy(), args = (Xdata,Gdata))
    return fit

# --- A function to crop down a data array --------
def crop(imgarray,xmin,xmax,ymin,ymax):
    """ Crop array 'imgarray' down to new limits defined by xmin,xmax,ymin,ymax.
 
    Beware of Python's array indexing convention!"""
    out = imgarray[xmin:xmax,ymin:ymax]
    return out

# --- 2D gaussian fit to 2D image data array ---------
def gfit2D( xdata,ydata,Gdata,p=array([1,5,5,10,10]) ):
    """ Fit 2D image data 'Gdata' to a 2D gaussian over the data range specified by 'xdata' 
    and 'ydata'.
        
        p is the array (5 values) of inital guesses:
        p[0] = amplitude
        p[1] = x width
        p[2] = y width
        p[3] = x centre
        p[4] = y centre """
    # 2D gaussian definition
    fitfunc = lambda p, X, Y: p[0]*exp( -(X-p[3])**2 /(2*(p[1]**2)) 
    -(Y-p[4])**2 / (2*(p[2]**2)) )
    # error function
    errfunc = lambda p, X, Y, val: fitfunc(p,X,Y) - val
    # perform the fit
    Xdata,Ydata = meshgrid(xdata,ydata) # create 2D 'x' and 'y' coord. arrays
    fit,success = optimize.leastsq(errfunc, p.copy(), args = (Xdata,Ydata,Gdata))
    # make plot of results - two cross-sections through centre of peak
    xrange = linspace(xdata.min(), xdata.max(),200)
    yrange = linspace(ydata.min(), ydata.max(),200)
    xcent = round(fit[3]); ycent = round(fit[4])    # rough array indices of cloud centres
    subplot(2,1,1)  # top left plot box (2 x 2 grid of axes 
    #fit,success = optimize.leastsq(errfunc, p.copy(), args = (Xdata,Ydata,Gdata) )
    plot(xrange,fitfunc(fit,xrange,yrange)[ycent,:],"k-",Xdata,Gdata[ycent,:],"b.") #plot data and fit
    subplot(2,1,2)
    plot(xrange,errfunc(fit,xrange,yrange,Gdata)[ycent,:],"b.-",Xdata,zeros((1,len(Xrange))),"b.") # plot residuals
    
#    ax = axes()
    #draw() fit,success = optimize.leastsq(errfunc, p.copy(), args = (Xdata,Ydata,Gdata) )
    show()
    return fit

def gaussian1D(x, params):
    """ Definition of 1D gaussian.  
        Call from withing an 'analysis' object with gaussian parameters 
        'params'.
        
        params[0] = amplitude
        params[1] = width**2 (sigma**2)
        params[2] = peak centre
        params[3] = constant background
            
        x = 1D dependent variable (e.g. indices of pixels)
    """
    out = params[0]*exp( -(x-params[2])**2 /(2*params[1]) ) + params[3]
    return out

def gaussian2D(X, Y, params):
    """ Definition of 2D gaussian.  
        Call from withing an 'analysis' object with gaussian parameters 'params'.
        
        params[0] = amplitude
        params[1] = x width squared (sigmaX**2)
        params[2] = y width squared (sigmaY**2)
        params[3] = x centre
        params[4] = y centre
        params[5] = constant background
            
        X,Y are 'meshgrid'-formatted 2D index arrays 
    """
    out = params[0]*exp( -(X-params[3])**2 /(2*params[1]) -(Y-params[4])**2 /(2*params[2]) ) + params[5]
    return out

def fit_prep(data_dic, matrix_view_object, splitting_object):
    """ Grab 1D data from the 2D image by ROI sum, or by cursor slice.
        Also compute and return pixel sum and mean backgroun level.
        Also use parameters in splitting_object to compute pixel sum in the
        two auxilliary boxes in 2D O.D. image.
    """
    mv = matrix_view_object
    sp = splitting_object
   
    # Integer-ize on-screen data limits into usable image array indices.
    # Parity between on-screen x indices and proper array indices thanks to
    # origin='lower' set in imshow() call in  update() method of MatrixView
    x_min = int(mv.x_min)
    x_max = int(mv.x_max) - 1   # Need the -1. mv.x,y_max always one too large.
    y_min = int(mv.y_min)
    y_max = int(mv.y_max) - 1
  
    # Set up dependent varaible (image data)
    if mv.data_subset == 'cursor':
        # Only fit to cursor data slices
        # Remember: we're working with potentially cropped data here.
        h_gdata = data_dic[mv.selected_image]\
                [mv.cursor_y_pos, x_min:x_max] 
        v_gdata = data_dic[mv.selected_image]\
                [y_min:y_max, mv.cursor_x_pos] 
    elif mv.data_subset == 'ROI sum':
        # Sum over data in current zoom range.
        data = data_dic[mv.selected_image][y_min:y_max, x_min:x_max]
        h_gdata = data.sum(axis=0)              
        v_gdata = data.sum(axis=1)              
    else:
        print "**** bad string passed for 'data_subset' in data_analysis.py "
    # Set up independent variable (coordinate vector)
    h_xdata = arange(x_min, x_max) 
    v_xdata = arange(y_min, y_max) 

    # Compute pixel sum, background level of total analysis region
    # Use a 3-pixel wide background perimeter of the image to define
    # background level.
    pix_sum, bg = get_pix_sum(3,\
            data_dic[mv.selected_image][y_min:y_max, x_min:x_max] )
     
    # Compute pixel sum, background level in left-hand analysis region
    # based on limits defined in litting_object. Use background 
    # computed from entire image.
    pix_sum_left = get_pix_sum2(bg,\
            data_dic[mv.selected_image]\
                [sp.left_vert_bot:sp.left_vert_top,\
                 sp.left_horiz_left:sp.left_horiz_right] )
    
    # Compute pixel sum, background level in left-hand analysis region
    # based on limits defined in splitting_object. Use background 
    # computed from entire image.
    pix_sum_right = get_pix_sum2(bg,\
            data_dic[mv.selected_image]\
                [sp.right_vert_bot:sp.right_vert_top,\
                 sp.right_horiz_left:sp.right_horiz_right] )

    return h_xdata, h_gdata, v_xdata, v_gdata, pix_sum, bg,\
            pix_sum_left, pix_sum_right

def fit_prep2(data_dic, matrix_view_object, splitting_object):
    """ Preparing image data from the left and right box (i.e clouds)
        Grab 1D data from the 2D image by ROI sum, or by cursor slice.
               
    """
    mv = matrix_view_object
    sp = splitting_object
   
##    # Integer-ize on-screen data limits into usable image array indices.
##    # Parity between on-screen x indices and proper array indices thanks to
##    # origin='lower' set in imshow() call in  update() method of MatrixView
##    x_min = int(mv.x_min)
##    x_max = int(mv.x_max) - 1   # Need the -1. mv.x,y_max always one too large.
##    y_min = int(mv.y_min)
##    y_max = int(mv.y_max) - 1
##  
##    # Set up dependent varaible (image data)
##    if mv.data_subset == 'cursor':
##        # Only fit to cursor data slices
##        # Remember: we're working with potentially cropped data here.
##        h_gdata = data_dic[mv.selected_image]\
##                [mv.cursor_y_pos, x_min:x_max] 
##        v_gdata = data_dic[mv.selected_image]\
##                [y_min:y_max, mv.cursor_x_pos] 
##    elif mv.data_subset == 'ROI sum':
##        # Sum over data in current zoom range.
##        data = data_dic[mv.selected_image][y_min:y_max, x_min:x_max]
##        h_gdata = data.sum(axis=0)              
##        v_gdata = data.sum(axis=1)              
##    else:
##        print "**** bad string passed for 'data_subset' in data_analysis.py "
##    # Set up independent variable (coordinate vector)
##    h_xdata = arange(x_min, x_max) 
##    v_xdata = arange(y_min, y_max) 
##
##    # Compute pixel sum, background level of total analysis region
##    # Use a 3-pixel wide background perimeter of the image to define
##    # background level.
##    pix_sum, bg = get_pix_sum(3,\
##            data_dic[mv.selected_image][y_min:y_max, x_min:x_max] )
##     
##    # Compute pixel sum, background level in left-hand analysis region
##    # based on limits defined in litting_object. Use background 
##    # computed from entire image.
##    pix_sum_left = get_pix_sum2(bg,\
##            data_dic[mv.selected_image]\
##                [sp.left_vert_bot:sp.left_vert_top,\
##                 sp.left_horiz_left:sp.left_horiz_right] )
##    
##    # Compute pixel sum, background level in left-hand analysis region
##    # based on limits defined in splitting_object. Use background 
##    # computed from entire image.
##    pix_sum_right = get_pix_sum2(bg,\
##            data_dic[mv.selected_image]\
##                [sp.right_vert_bot:sp.right_vert_top,\
##                 sp.right_horiz_left:sp.right_horiz_right] )

  

    #Defining the limits of the boxes
   
    L_xmin = sp.left_horiz_left
    L_xmax = sp.left_horiz_right
    R_xmin = sp.right_horiz_left
    R_xmax = sp.right_horiz_right

    L_ymin = sp.left_vert_bot
    L_ymax = sp.left_vert_top
    R_ymin = sp.right_vert_bot
    R_ymax = sp.right_vert_top

    
    #The "x" independent data used for 1D fits. We do 1D fits on the vertical
    #and horizontal data for both left and right boxes "i.e clouds". We may need to subtract / add +1...
    #horizontal coordinate data arrays for left and right boxes

    L_h_xdata = arange(sp.left_horiz_left, sp.left_horiz_right)
    R_h_xdata = arange(sp.right_horiz_left, sp.right_horiz_right)

    #vertical coordinate data arrays for left and right boxes
    L_v_xdata = arange(sp.left_vert_bot, sp.left_vert_top)
    R_v_xdata = arange(sp.right_vert_bot, sp.right_vert_top)


    #Getting the image data arrays for each horizontal and L/R
    # Each point integrates the transverse direction.

    L_data = data_dic[mv.selected_image][L_ymin:L_ymax, L_xmin:L_xmax]
    R_data = data_dic[mv.selected_image][R_ymin:R_ymax, R_xmin:R_xmax]
    
    L_h_gdata = L_data.sum(axis=0)              
    L_v_gdata = L_data.sum(axis=1)

    R_h_gdata = R_data.sum(axis=0)              
    R_v_gdata = R_data.sum(axis=1)

         

    return L_h_xdata, R_h_xdata, L_v_xdata, R_v_xdata, L_h_gdata,\
           L_v_gdata, R_h_gdata, R_v_gdata


def get_pix_sum(N, img):
    """ Compute pixel sum of optical density in given data array 'img'.
        Subtract off background first by computing average OD in a
	perimeter N pixels thick at the image edges.
    """
    # The perimeter bits of data
    left = (img[:,0:N]).mean()
    right = (img[:,-N:]).mean()
    top = (img[0:N,N:-N]).mean()
    bot = (img[-N:,N:-N]).mean()
    # mean background
    bg = (left + right + top + bot)/4
    # Return pixel sum and mean background
    return (img-bg).sum(), bg

def get_pix_sum2(bg, img):
    """ Compute pixel sum of optical density in given data array 'img'.
        Subtract off background specified by 'bg'.  Use this for
        left / right pixel sums.
    """ 
    # Return pixel sum and mean background
    return (img-bg).sum()
    
def get_peak_OD(data_dic, matrix_view_object, h_fit, v_fit, bg):
    """ Access the 2D data array using peak coordinates determined by
        fits.  Return peak O.D. value (corrected for background!)
    """
    mv = matrix_view_object
    
    # Integer-ize on-screen data limits into usable image array indices.
    # (Copied from fit_prep method)
    x_min = int(mv.x_min)
    x_max = int(mv.x_max) - 1   # Need the -1. mv.x,y_max always one too large.
    y_min = int(mv.y_min)
    y_max = int(mv.y_max) - 1
    
    # The working 2D data set.
    data = data_dic[mv.selected_image][y_min:y_max, x_min:x_max]
    # Access data array at fit peak centres.  Subtract off mean bkgnd.
    peak_OD = data[v_fit[2], h_fit[2]] - bg
    return peak_OD

def seed_gaussian1D(h_xdata, h_gdata, v_xdata, v_gdata, matrix_view_object):
    """ Cleverly guestimate seed parameters for 1D gaussian fit.
    """
    mv = matrix_view_object

    #######################################################################
    # To better estimate fit parameters, smooth data first. (4-pix. window)
    #######################################################################
    ind1h = arange(0, h_gdata.shape[0])
    h_gsmooth = ( h_gdata[ind1h] + h_gdata[ind1h-1] + 
                  h_gdata[ind1h-2] + h_gdata[ind1h-3] ) /4

    ind1v = arange(0, v_gdata.shape[0])
    v_gsmooth = ( v_gdata[ind1v] + v_gdata[ind1v-1] + 
                  v_gdata[ind1v-2] + v_gdata[ind1v-3] ) /4
    # Estimate background level from outer 12 pixels.
    # Subtract two since smoothed peak should be shifted over by 1.5 pix.
    h_bg = (h_gsmooth[0:6] + h_gsmooth[-6:]).sum()/12 - 1.5
    v_bg = (v_gsmooth[0:6] + v_gsmooth[-6:]).sum()/12 - 1.5
    # Find ind. of those points which are greater than half the peak height.
    # Hopefully h,v_smooth.max() roughly corresponds to the real peak.
    # Subtract off backgrounds first.
    h_gs = h_gsmooth - h_bg
    v_gs = v_gsmooth - v_bg
    ind2h = find( h_gs > h_gs.max()/2 )
    ind2v = find( v_gs > v_gs.max()/2 ) 
    #print "ind2h size = %d " % ind2h.shape[0]
    #print "ind2v size = %d " % ind2v.shape[0]
    # Estimate peak centres using barycenter (C.O.Mass) of smoothed data.
    # h_c and v_c will refer to real data indices, not the relative 'ind2'.
    try:
        h_c = (h_xdata[ind2h]*h_gs[ind2h]).sum() /\
                float(h_gs[ind2h].sum())
    except Exception, inst:
        # Set centre to middle index in case of out-of-bounds error.
        h_c = h_xdata[h_xdata.shape[0]/2]
        print "******* peak indices out of bounds (horiz.) ********"
    try:
        v_c = (v_xdata[ind2v]*v_gs[ind2v]).sum() /\
                float(v_gs[ind2v].sum())
    except Exception, inst:
        # Set centre to middle index in case of out-of-bounds error.
        v_c = v_xdata[v_xdata.shape[0]/2]
        print "******* peak indices out of bounds (vert.) ********"
    #print ('x_min', x_min, 'x_max',x_max, 'y_min', y_min, 'y_max', y_max) 

    # Estimate amplitudes - max value among 'ind2' 
    h_ampd = h_gs[ind2h].max()
    v_ampd = v_gs[ind2v].max()

    # Estimate widths-squared - second moment of peak of smoothed.
    # Throw in a factor of 2**2 for good measure.
    h_sigma2 = ((h_xdata[ind2h]-h_c)**2 * abs(h_gs[ind2h])).sum() /\
            float(abs(h_gs[ind2h]).sum())
    v_sigma2 = ((v_xdata[ind2v]-v_c)**2 * abs(v_gs[ind2v])).sum() /\
            float(abs(v_gs[ind2v]).sum())
    
    return (array([h_ampd, h_sigma2, h_c, h_bg]), 
                array([v_ampd, v_sigma2, v_c, v_bg]))

def crunch_split_params(experiment_object, camera_object, pix_sum_left,\
            pix_sum_right):
    """ Compute atom number in left/right regions given pixel sum info,
        as well as fraction in left well.
    """
    # Pointers to the experiment and camera objects.
    expt = experiment_object
    cam = camera_object
    # Define physical constants
    # Note: X = horizontal, Y = vertical, 
    #       Z = integrated (i.e. into/out of image)
    ###################################################
    hbar = 1.05459e-34              # [J/s]
    kb = 1.3807e-23                 # [J/K]
    tof = 0.001* expt.tof           # time of flight [s]
    umX = 2.72e-6 * cam.binning_X   # Effective pixel size, horiz. and vert.
    umY = 2.72e-6 * cam.binning_Y   # 13 [um/pix] / 4 = 3.25 [um/pix] 
    
    # Define atomic mass [kg], natural linewidth [Hz] and wavelength [m].
    if expt.atoms == 'Rb':
        m = 1.44e-25
        gamma = 5.98e6 # Real frequency units [Hz] 
        wlen = 780.027e-9
    elif expt.atoms == 'K':
        m = 6.67e-26
        gamma = 6e6 # FIXME: get the correct value here
        wlen = 767.7017e-9
    else:
        print "**** bad string passed for experiment.atom ****"
    # Imaging beam detuning in units of half-natural-linewidhts
    det = expt.detuning*1e6 / (gamma/2)
    # Absorption cross-section
    sigma_absn = ((3*(wlen)**2)/(2*pi)) / (1 + det**2)
    # Compute the pixel sum atom numbers and left fraction 
    N_ps_left = pix_sum_left*umX*umY / sigma_absn
    N_ps_right = pix_sum_right*umX*umY / sigma_absn
    p_L = N_ps_left / (N_ps_left + N_ps_right)
    
    return N_ps_left, N_ps_right, p_L

def crunch_params(experiment_object, camera_object, matrix_view_object,\
            h_fit, v_fit, L_h_fit, L_v_fit, R_h_fit, R_v_fit,\
            pix_sum,  pix_sum_left,\
            pix_sum_right, param_list, fit_type="g"):
    """ Function which converts results of the data fit into relevant
        physical parameters to be displayed and (possibly) saved.
        Gaussian fit (fit_type == "g"):
            h_fit[0] = amplitude
            h_fit[1] = width squared (sigma**2)
            h_fit[2] = centre
            h_fit[3] = const. background
        Same for v_fit. See "gfit1D()".
        'pix_sum' is the 2D pixel sum
    """
    # Pointers to the experiment, camera and matrix_view objects.
    expt = experiment_object
    cam = camera_object
    mv = matrix_view_object
    # Define physical constants
    # Note: X = horizontal, Y = vertical, 
    #       Z = integrated (i.e. into/out of image)
    ###################################################
    hbar = 1.05459e-34              # [J/s]
    kb = 1.3807e-23                 # [J/K]
    tof = 0.001* expt.tof           # time of flight [s]
    umX = 2.72e-6 * cam.binning_X   # Effective pixel size, horiz. and vert.
    umY = 2.72e-6 * cam.binning_Y   # 13 [um/pix] / 4 = 3.25 [um/pix] 
    # Define trap frequencies and anisotropy parameter 'L' 
    if expt.view_point == 'axial':
        # axial imaging
        omegaX = 2*pi*expt.trap_freq_ax
        omegaY = 2*pi*expt.trap_freq_rad
        omegaZ = 2*pi*expt.trap_freq_rad
        L = sqrt(omegaX**2 / omegaZ**2) # CHECKME: not sure about this!!
        omega_bar = (omegaX*omegaY*omegaZ)**(.3333333) #Added by Nathan, May 2012
    elif expt.view_point == 'radial':
        # radial imaging
        omegaX = 2*pi*expt.trap_freq_rad
        omegaY = 2*pi*expt.trap_freq_rad
        omegaZ = 2*pi*expt.trap_freq_ax
        L = sqrt(omegaZ**2 / omegaX**2)
        omega_bar = (omegaX*omegaY*omegaZ)**(.3333333) #Added by Nathan, May 2012
    else:
        print "**** bad string passed for experiment.view_point ****"
    # Define atomic mass [kg], natural linewidth [Hz] and wavelength [m].
    if expt.atoms == 'Rb':
        m = 1.44e-25
        gamma = 5.98e6 # Real frequency units [Hz] 
        wlen = 780.027e-9
    elif expt.atoms == 'K':
        m = 6.67e-26
        gamma = 6e6 # FIXME: get the correct value here
        wlen = 767.7017e-9
    else:
        print "**** bad string passed for experiment.atom ****"
    # Imaging beam detuning in units of half-natural-linewidhts
    det = expt.detuning*1e6 / (gamma/2)
    # Absorption cross-section
    sigma_absn = ((3*(wlen)**2)/(2*pi)) / (1 + det**2)
    # Compute the physical parameters
    ###################################################
    # Build disctionary of phys. parameters. Initialize to 0.
    phys_params = {}
    for key in param_list:
        # Exclude the "fit" entry - it's not a real parameter.
        if key != "fit":
            phys_params[key] = 0
    # Gaussian fit analysis
    if fit_type == 'g':
        # Gaussian widths after TOF [pixels]
        phys_params['sigmaX'] = sqrt(h_fit[1])
        phys_params['sigmaY'] = sqrt(v_fit[1])
        # Compute atom number for 'ROI sum' or from 'cursor' fits:
        if mv.data_subset == 'ROI sum':
            # Atom number from 1D fit to integrated O.D. data 
            # USE "ROI SUM" FOR ACCURATE ATOM NUMBER (fit)
            # NB. Pixel area is 1 x sigma_h pixels.  Thus, multiply by
            # umX*umY to convert to area in metres-squared.
            phys_params['N_fit'] = \
                h_fit[0]*sqrt(2*pi)*phys_params['sigmaX']*umX*umY / sigma_absn
        elif mv.data_subset == 'cursor':
            phys_params['N_fit'] = \
                2*pi*h_fit[0] *(phys_params['sigmaX'] *umX)**2 / sigma_absn
        else:
            print "**** bad string passed for mv.data_subset ***"
            phys_params['N_fit'] = -1

        #Pixel sum (O.D. sum) atom number
        phys_params['N_pix_sum'] = pix_sum*umX*umY / sigma_absn
        # Temperatures [Kelvin]
        phys_params['Tx'] = \
            ( m*omegaX**2 *(h_fit[1])*(umX**2) )/( kb*(1+omegaX**2 * tof**2) )
        phys_params['Ty'] = \
            ( m*omegaY**2 *(v_fit[1])*(umY**2) )/( kb*(1+omegaY**2 * tof**2) )

        #PEAK CENTER!

        
        # Peak centres [pixels]
        phys_params['x_centre'] = h_fit[2]
        phys_params['y_centre'] = v_fit[2]
        phys_params['L_y'] = L_v_fit[2]
        phys_params['R_y'] = R_v_fit[2]
        phys_params['L_x'] = L_h_fit[2]
        phys_params['R_x'] = R_h_fit[2]
        phys_params['delta_y'] = L_v_fit[2] - R_v_fit[2]
        phys_params['delta_x'] = L_h_fit[2] - R_h_fit[2]
        
               
        # Fermi temperature
        Tf = hbar*omega_bar*( 6*(phys_params['N_fit']) )**(0.333333) / kb
        phys_params['TxoTf'] = phys_params['Tx'] / Tf
        # Inferred in-trap gaussian widths (TOF = 0)
        sigmaX0 = sqrt(kb*phys_params['Tx'] /(m*omegaX**2))
        sigmaY0 = sqrt(kb*phys_params['Ty'] /(m*omegaY**2))
        Tz = sqrt(phys_params['Tx'] * phys_params['Ty'])
        sigmaZ0 = sqrt(kb*Tz /(m*omegaZ**2))
        # Peak phase space density in trap
        n_0 = phys_params['N_fit'] / (sigmaX0*sigmaY0*sigmaZ0 *(2*pi)**1.5)
        ldb = sqrt( (2*pi*hbar**2) / (m*kb*Tz) ) 
        phys_params['psd'] = n_0*(ldb**3)
        # Compute Tc
        omega_bar = (omegaX*omegaY*omegaZ)**(0.333333)
        Tc = (0.94*hbar*omega_bar*(phys_params['N_fit'])**(0.333333) )/kb
        # Print Tc and n_0 to "standard output"
        print "Tc = %3.3e K" % Tc
        print "peak density = %3.3e m-3" % n_0
        phys_params['N_L'] = pix_sum_left*umX*umY / sigma_absn
        phys_params['N_R'] = pix_sum_right*umX*umY / sigma_absn
        phys_params['p_L'] = \
            phys_params['N_L'] / (phys_params['N_L'] + phys_params['N_R'])
        phys_params['z'] =  \
           (phys_params['N_R'] - phys_params['N_L']) / (phys_params['N_L'] + phys_params['N_R'])  
    elif fit_type == 'f':
        print "do some fermi fit analysis here"
    else: print "**** bad string for fit_type *****"
    # Check for values which are too large or too small for matlab
    #for key, value in phys_params.iteritems():
    #if abs( log( abs(value) ) ) > 18:    # the exponent
    #        phys_params[key] = -1
    #
    return phys_params
    
def pretty_print(item):
    """ Returns a string representation of an item.
    """
    if isinstance(item, str):
        return item
    elif isinstance(item, float):
        if item > 1e3 or item < -1e3 or (item < 1e-2 and item > -1e-2):
            return "%.3e" % item
        else:
            return "%.3f" % item
    elif isinstance(item, int):
        if item > 9999 or item < - 9999:
            return "%.0e" % item
        else:
            return "%d" % item
    else:
        return item.__repr__()

def matlab_format(item):
    """ Returns decimal notation float. of item.
        Format values in this way so that they are loadable in Matlab
        after being exported with 'io.savemat' method. """

    temp_str = "%.10f" % item
    print temp_str
    return float(temp_str)

def dic_to_string(values_dic, units_dic):
    """ Returns a nicely formatted string from two dictionaries. """
    print_key = lambda key : "%s: %s %s" % (pretty_print(key),\
            pretty_print(values_dic[key]), pretty_print(units_dic[key]) )
    return "\n".join((print_key(key) for key in values_dic))

def display_params(dic, ordered_units_list):
    """ Format text output for physical parameters.
        Parameters are stored in dictionary 'dic', with keys given by
        'param_list'
    """
    print_key = lambda key: "%s: %s %s" % ( pretty_print(key),\
        pretty_print(dic[key]), pretty_print(ordered_units_list[key]) )
    # string type output
    return "\n".join((print_key(key) for key in dic))

def make_OD(a_image, a_bg, l_image, l_bg):
    """ Construct a divided optical density array using the given input
        image arrays.  'a', 'l' = 'atom' and 'laser' images (i.e. no atoms).
        'a_bg' and 'l_bg' are background levels to be subtracted before making
        divided image.
    """
    thresh = 30
    a = a_image - a_bg
    l = l_image - l_bg

    # Fill in zero pixels with ones to avoid div-by-zero problems.
    absn = (a + 1*(a==0)) / (l + 1*(l==0))

    return -log(absn) 

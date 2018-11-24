import os
import sys
import numpy as np
import matplotlib.pyplot as pp
import warnings
from datetime import datetime

try:
    from PIL import Image
    PIL_imported = True
except ImportError:
    PIL_imported = False
    print "Python module 'PIL' not available."
    pass



class EOS1_Img:
    """Python class for analyzing spectra images taken with the EOS 1 open-source spectrometer.
    
    Jiansheng Feng  ( jfeng@uakron.edu  www.erieopen.tech )
    
    To initialize, need input: 
        - file name of the image
        - file path (from the Python folder)
        - whether to trim edge #default False#
        
    Functions:
        - show_RGB
        - show_hmap: requires running cal_heatmap first
        - cal_heatmap
        - find_ref: requires running cal_heatmap first
        - check_align: requires running find_ref (and therefore cal_heatmap) first
        - norm_sam: requires running find_ref (and therefore cal_heatmap) first
        - cal_I: requires running norm_sam (and therefore find_ref and cal_heatmap) first
        - test_N: runs cal_heatmap, find_ref, norm_sam, cal_I; requires nitrate_calibration.csv
    """
    
    
    def __init__( self, filename, filepath, trim_edge=False ):
        self.filename = filename
        self.filepath = filepath
        
        self.xi = pp.imread( self.filepath+self.filename )
        if self.xi.ndim != 3:
            raise ValueError( "Image array dimension incorrect. Expect an RGB image." ) 
        else:
            pass
        if self.xi.dtype != np.uint8:
            raise ValueError( "Image array datatype incorrect. Expect an 8-bit image." ) 
        else:
            pass
        
        self.h, self.w, self.d = self.xi.shape 
        if self.h < self.w:
            warnings.warn( "Image appears to be landscape." )
            proceed = raw_input( "Continue? (y/N): " )
            if proceed=='y' or proceed=='Y':
                pass
            else:
                raise RuntimeError( "Program terminated by user." )
        else:
            pass
        
        if trim_edge == True:
            self.xi = self.xi[self.h/4:self.h*3/4, self.w/4:self.w*3/4, :]
            self.show_RGB( fig_dpi=100 )
            proceed = raw_input( "Continue with this image? (y/N): " )
            if proceed=='y' or proceed=='Y':
                pass
            else:
                raise RuntimeError( "Restart and set 'trim_edge' to False." )
        else:
            pass
        
    
    def show_RGB( self, fig_dpi=200 ):
        pp.figure( dpi=fig_dpi )
        pp.style.use( "seaborn-dark" )
        pp.imshow( self.xi )
        
        
    def show_hmap( self, fig_dpi=200, x_max=510 ):
        pp.figure( dpi=fig_dpi )
        pp.style.use( "seaborn-dark" )
        pp.imshow( self.xh, cmap="gray", vmin=0, vmax=x_max )
        pp.colorbar()
        
    
    def cal_heatmap( self ):
        self.xf = self.xi.astype( float )
        self.xh = abs( self.xf[:,:,0] - self.xf[:,:,1] )
        self.xh += abs( self.xf[:,:,0] - self.xf[:,:,2] )
        self.xh += abs( self.xf[:,:,1] - self.xf[:,:,2] )
        
        
    def find_ref( self, n=0.25 ):
        self.n = float( n )
        if self.n<0.1 or self.n>0.9:
            self.n = 0.25
        else:
            pass
        
        self.x0 = self.xh.mean( axis=0 )
        self.x0thres = np.argwhere( self.x0 > self.x0.max()*n ).flatten()
        self.x0diff  = self.x0thres[1:] - self.x0thres[:-1]
        self.x0gap   = np.where( self.x0diff > 2. )[0].flatten()
        if len( self.x0gap )==0:
            self.l_edge, self.r_edge = self.x0thres[0], self.x0thres[-1]
        else:
            self.d_to_center = []
        for i in self.x0gap:
            self.d_to_center.append( abs( self.w/2. - self.x0thres[i:i+2].mean() ) )
        self.d_min = np.argmin( self.d_to_center )
        if self.d_min==0:
            self.l_edge, self.r_edge = self.x0thres[0], self.x0thres[ self.x0gap[0] ]
        else:
            self.l_edge = self.x0thres[ self.x0gap[self.d_min-1]+1 ]
            self.r_edge = self.x0thres[ self.x0gap[self.d_min] ]
            
        self.x_ref_band = self.xh[ :, self.l_edge:self.r_edge+1 ]
        self.x1 = self.x_ref_band.mean( axis=1 )
        self.x1thres = np.argwhere( self.x1 > self.x1.max()*n ).flatten()
        self.t_edge, self.b_edge = self.x1thres[0], self.x1thres[-1]
        
        self.ref_wid = self.r_edge - self.l_edge
        self.half_hgt = int( (self.b_edge - self.t_edge)/2. )
        
    
    def check_align( self ):
        xh_edge = self.xh[self.t_edge:self.b_edge+1, 
                          self.r_edge-self.ref_wid/2:self.r_edge+self.ref_wid/2]
        xhslpmx = []
        for row in range(len(xh_edge)):
            xhslope = xh_edge[row, 1:]-xh_edge[row, :-1]
            xhslpmx.append( np.argmin( xhslope ) )
        xhslpmx = np.array( xhslpmx )
        self.fedge = np.polyfit( np.arange(len(xhslpmx)), xhslpmx, 1 )[0]
        if abs(self.fedge) > 0.1:
            print "\nImage tilt ~ "+str( np.arctan(self.fedge)*180./3.1416 )+" outside permitted range."
            rotate = raw_input( "Let program rotate image? (Y/n):" )
            if rotate == 'N' or rotate == 'n': 
                print ( "If continue with this image, it could result in significant measurement error." )
                proceed = raw_input( "Continue anyway? (y/N): " )
                if proceed=='y' or proceed=='Y':
                    pass
                else:
                    raise RuntimeError( "Please rotate image manually and restart." )
            elif PIL_imported == False:
                raise RuntimeError( "PIL module not imported. Please rotate image manually and restart." )
            else:
                img = Image.open( self.filepath+self.filename )
                file_ID, file_ext = self.filename.split('.')
                self.filename_rot = file_ID+"_rot."+file_ext
                img_rot = img.rotate( -1.*np.arctan(self.fedge)*180./3.1416 )
                img_rot.save( self.filepath+self.filename_rot )
                print ( "Rotated image "+self.filename_rot+" generated and saved." )
                self.xi = pp.imread( self.filepath+self.filename_rot )
                self.h, self.w, self.d = self.xi.shape
                self.cal_heatmap()
                self.find_ref() 
        else:
            print "\nImage tilt ~ "+str( np.arctan(self.fedge)*180./3.1416 )
            print "Alignment check passed."
        
        
    def norm_sam( self, bpeak_chl='r', trim_margin=True, gapcal='p' ):
        if trim_margin == True:
            self.mrg = int( self.ref_wid/10. )
        else:
            self.mrg = 0
        
        self.x_ref = self.xi[ self.t_edge:self.b_edge, 
                              self.l_edge+self.mrg:self.r_edge-self.mrg, : ]
        self.y_ref = self.x_ref.mean( axis=1 )
    
        if bpeak_chl == 'r':
            self.peak_r = self.y_ref[:self.half_hgt,0].argmax()
            self.peak_b = self.y_ref[self.half_hgt:,0].argmax()+self.half_hgt
        else:
            self.peak_rgb = self.y_ref.argmax( axis=0 )
            self.peak_r, self.peak_b = self.peak_rgb[[0,2]]

        if gapcal == 'w':
            self.gap = int( self.ref_wid*0.901 )
        else:
            self.gap = int( ( self.peak_b-self.peak_r )*0.368 )
    
        self.x_sam = self.xi[ self.t_edge:self.b_edge, 
                              self.r_edge+self.gap+self.mrg:self.r_edge+self.gap+self.ref_wid-self.mrg, : ]
        self.y_sam = self.x_sam.mean( axis=1 )
        self.max_rgb = self.y_ref.max( axis=0 )

        self.peak_px = np.array([self.peak_r, self.peak_b]).flatten()
        self.peak_nm = np.array([610.65, 449.1])
        self.f = np.polyfit( self.peak_px, self.peak_nm, 1 )
        self.wl = np.arange(self.b_edge-self.t_edge)*self.f[0]+self.f[1]
    
        self.y_sam_norm_r = self.y_sam[:, 0]/self.max_rgb[0]
        self.y_sam_norm_g = self.y_sam[:, 1]/self.max_rgb[1]
        self.y_sam_norm_b = self.y_sam[:, 2]/self.max_rgb[2]
        self.y_sam_norm = np.dstack((self.y_sam_norm_r, self.y_sam_norm_g, self.y_sam_norm_b))[0]


    def cal_I( self, chl='g', wl_low=525., wl_high=535. ):  
        if chl=='r' or chl=='R':
            self.sam_arr = self.y_sam_norm_r
        elif chl=='g' or chl=='G':
            self.sam_arr = self.y_sam_norm_g
        elif chl=='b' or chl=='B':
            self.sam_arr = self.y_sam_norm_b
        else:
            raise ValueError( "Color channel should be 'r', 'g', or 'b'." )
        
        arg_low = np.where( self.wl < wl_high )[0][0]
        arg_high = np.where( self.wl > wl_low )[0][-1]
        I_sum = self.sam_arr[arg_low:arg_high+1].sum()
        I_ave = I_sum/(arg_high-arg_low+1)
        return I_ave
        
    
    def test_N( self, wlc=530., wlhw=5., cali="nitrate_calibration.csv" ):
        wl_l = wlc-wlhw
        wl_h = wlc+wlhw
        
        self.cal_heatmap()
        self.find_ref()
        self.check_align()
        self.norm_sam()
        I_ave = self.cal_I( wl_low=wl_l, wl_high=wl_h )        
        
        if os.path.isfile(cali) == True:
            f = open( cali )
            ls = f.readlines()
            f.close()
            cali_date = ls[0].strip().split(',')[1]
            cali_time = ls[1].strip().split(',')[1] 
            k_N = float( ls[2].strip().split(',')[1] )
            b_N = float( ls[3].strip().split(',')[1] )
            print "\nUsing calibration record from "+cali_date+" "+cali_time
        else:
            print "Calibration record not found."
            input_kb = raw_input( "\n\nManually input k and b values? (y/N): " )
            if input_kb == 'y' or input_kb == 'Y':
                k_N = float( raw_input("Please specify k_N: ") )
                b_N = float( raw_input("Please specify b_N: ") )
            else:
                k_N = -7.8279
                b_N = -0.14917
        
        lgI = np.log10(I_ave)
        nc = lgI*k_N + b_N
        return nc
        
    

def cali_N( img_arr, nc_arr, f_path, wl=530., fig_out=True ):
    if len(img_arr) != len(nc_arr):
        raise ValueError( "img_arr and nc_arr should have the same length." )
    else:
        pass
    nc_arr = np.array( nc_arr )
    
    I_arr = []
    for img in img_arr:
        eos = EOS1_Img( img, f_path )
        eos.cal_heatmap()
        eos.find_ref() 
        eos.check_align() 
        eos.norm_sam()
        I_arr.append( eos.cal_I('g', wl-5., wl+5.) )
    I_arr = np.array( I_arr )
    lgI = np.log10(I_arr)    
    k, b = np.polyfit( lgI, nc_arr, 1 )
    
    cali_rec = "nitrate_calibration.csv"
    if os.path.isfile( cali_rec ) == True:
        print "Calibration record '"+cali_rec+"' already exist."
        proceed = raw_input("Overwrite this record? (Y/n): ")
        if proceed == 'n' or proceed == 'N':
            cali_rec = raw_input("File name of new record (including file extension): ")
        else:
            pass
    else:
        pass
    
    print "Writing to calibration record: "+cali_rec+'.'
    f = open( cali_rec, 'w' )
    day_str, time_str = str( datetime.now() ).split()
    time_str = time_str.split('.')[0]
    f.write("date,"+day_str+'\n' )
    f.write("time,"+time_str+'\n' )
    f.write("k,"+str(round(k, 5))+'\n' )
    f.write("b,"+str(round(b, 5))+'\n' )
    f.close()
    print "Done writing to calibration record."
    
    if fig_out == True:
        Ab_arr = (-1.)*lgI
        kf, bf = np.polyfit( nc_arr, Ab_arr, 1 )

        pp.style.use( "seaborn-darkgrid" )
        pp.figure( dpi=150 )
        pp.plot( nc_arr, Ab_arr, 'k.', label="Calibration Data" )
        pp.plot( nc_arr, nc_arr*kf+bf, 'k-', label="Linear Fit" )
        pp.xlabel( "Nitrate Concentration (mg/L)", size=12)
        pp.ylabel( "Absorbance ("+str(wl-5)+"nm $-$ "+str(wl+5)+"nm)", size=12 )
        pp.legend( loc="upper left" )
        pp.show()
    else:
        pass
    
        
        
if __name__ == "__main__":
    print "###############"
    print "This script is written for Python v2.7.15"
    print "You are using Python "+str( sys.version )
    print "###############"
    print '\n'
    print "Please make sure: "
    print "-1- The picture is portrait (i.e., height > width)."
    print "-2- The reference spectrum is on the left-hand-side."
    print "-3- The images are 8-bit RGB images (usually .jpg or .tif, but not .png)."
    print '\n'
    
    update_cali = raw_input( "Generate or update calibration record? (y/N): " )
    if update_cali=='y' or update_cali=='Y':
        img_path = raw_input( "Please specify the path to the images: " ) 
        img_list_str = raw_input( "Please list images, separated by commas: " ) 
        nc_list_str = raw_input( "Please list concentrations, separated by commas: " )
        wavelength = float( raw_input( "Please specify wavelength (nm): " ) )
        
        img_list, nc_list = [], []
        for img_str in img_list_str.strip().split(','):
            img_list.append( img_str.strip() )
        for nc_str in nc_list_str.strip().split(','):
            nc_list.append( float( nc_str.strip() ) )
        
        cali_N( img_list, nc_list, img_path, wavelength )
    else:
        pass
    
    meas_N = raw_input( "Measure nitrate? (y/N): " )
    if meas_N=='y' or meas_N=='Y':
        file_path = raw_input( "Please input image file path: " ) 
        file_name = raw_input( "Please input image file name (include extension): " ) 
        wavelength = float( raw_input( "Please input wavelength (nm): " ) )
        eos = EOS1_Img( file_name, file_path )
        nc = eos.test_N( wlc=wavelength )
        print "Nitrate Concentration: "+str(round(nc, 2))+" mg/L" 
    else:
        pass
    
    

#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# created by Jan Elias
# jan.elias@vut.cz
# Brno University of Technology,
# 2020

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import norm, gamma, weibull_min, lognorm, uniform, triang, truncnorm
#from scipy import sparse
from scipy.spatial import Voronoi, Delaunay
from scipy.special import gamma as gammafunc
import shutil
from freecad.woodWorkbench.tools.grafted import GD, findGD
from freecad.woodWorkbench.tools.custom_dist import pointDist
from scipy.stats.stats import pearsonr, spearmanr
from scipy.interpolate import interp1d

############################################################################################################    
############################################################################################################  
def crateDistributionObject(disttype,distparams):
    dist = []
    dist_stat = []
    if disttype == "Gaussian": 
        #mean + standard deviation
        dist = norm(distparams[0], distparams[1]) 
        dist_stat = [distparams[0],distparams[1]]
    if disttype == "TruncatedGaussian": 
        #mean, standard deviation, cuttoff        
        dist = truncnorm((distparams[2]-distparams[0])/distparams[1], float('inf'), loc=distparams[0], scale=distparams[1]) 
        mean, variance = dist.stats()
        dist_stat = [mean, variance**0.5]
    elif disttype == "Triangular": 
        #a, b, c (peak position)
        dist = triang((distparams[2]-distparams[0])/(distparams[1]-distparams[0]), distparams[0], distparams[1]-distparams[0])        
        mean, variance = dist.stats()
        dist_stat = [mean, variance**0.5]
    elif disttype == "Uniform": 
        #a + b
        dist = uniform(distparams[0], distparams[1]-distparams[0])        
        dist_stat = [0.5*(distparams[0]+distparams[1]),1./np.sqrt(12.)*(distparams[1]-distparams[0])]
    elif disttype == "Gamma": 
        #mean + standard deviation
        shape = (distparams[0]/distparams[1])**2
        scale = (distparams[1]**2)/distparams[0]
        dist = gamma(shape, 0., scale)
        dist_stat = [distparams[0],distparams[1]]
    elif disttype == "Weibull": 
        #mean + shape
        mean = distparams[0]
        shape = distparams[1]
        scale = mean/gammafunc(1.+1./shape)
        dist = weibull_min(shape, 0., scale)
        mean, variance = dist.stats()
        dist_stat = [mean, variance**0.5]
    elif disttype == "Lognormal": 
        #mean + standard deviation
        mean = distparams[0]
        std = distparams[1]
        s = (np.log(1.+(std/mean)**2))**0.5
        scale = np.exp( np.log((mean**2) / ((mean**2 + std**2)**0.5) ) )
        dist = lognorm(s, 0., scale)
        dist_stat = [mean,std]
    elif disttype == "Grafted": 
        dist =  findGD(mean = distparams[0], std = distparams[1], m = distparams[2], Pgr =distparams[3]) 
        dist_stat = [distparams[0],distparams[1]]
    elif disttype == "file":
        filename=distparams[0]
        d = np.loadtxt(filename)
        dist = pointDist(d)
        dist_stat = [dist.mean,dist.std]
    else: raise ValueError('Unknown distribution type: "%s"'%disttype)
    return dist, dist_stat

############################################################################################################    
############################################################################################################    

class RandomField():    
    def __init__(self, dimension = 3, readFromFolder = "none", dist_types=["Gaussian"], CC = np.array([[1]]), dist_params=[[0,1]], \
                 name="FIELD",out="OUTDIR",corr_l = 0.5, corr_f='exponential_square',x_range=[0.,1],y_range=[0.,2],z_range=[0.,0.5], sampling_type = "Sobol", filesavetype="binary"):

        if readFromFolder == "none":
            self.dim = dimension            

            self.cross_correlation_matrix =np.array(CC)
            self.Nvar = self.cross_correlation_matrix.shape[0]

            if( (not self.Nvar == len(dist_types)) or (not self.Nvar == len(dist_params)) ):
                raise ValueError("Error: number of random field in 'dist_types', 'dist_params' and 'CC' does not correspond")

            self.name = name
        
            self.filesavetype = filesavetype

            self.dist_types = dist_types
            self.dist_params = dist_params

            self.corr_l = corr_l #correlation length
            self.corr_f = "square_exponential" #correlation function

            self.x_range = np.array(x_range)
            self.y_range = np.array(y_range)
            if self.dim<2: self.y_range=np.array([0,0])
            self.z_range = np.array(z_range)
            if self.dim<3: self.z_range=np.array([0,0])

            self.sampling_type = sampling_type# MC, LHS, Freet or Sobol
            self.folder = Path(out + '/' + name + '/' + 'randomfield')
            self.saveReport()
        else: 
            self.folder = readFromFolder
            self.loadReport()
            
        self.CreateDistributionObjects()
        self.fraction_trace_considered = 0.99
        self.corr_l2 = self.corr_l**2
        self.grid_spacing = corr_l/3.;#for computational reason might be /3
        self.gap_mult = 1. #multipler of corr. length around the box to remove boundary bias
        self.x_range_gap = np.array([x_range[0]-self.gap_mult*corr_l,x_range[1]+self.gap_mult*corr_l])
        self.y_range_gap = np.array([y_range[0]-self.gap_mult*corr_l,y_range[1]+self.gap_mult*corr_l])
        if self.dim<2: self.y_range_gap=np.array([0,0])
        self.z_range_gap = np.array([z_range[0]-self.gap_mult*corr_l,z_range[1]+self.gap_mult*corr_l])
        if self.dim<2: self.z_range_gap=np.array([0,0])
        self.covmat_cutoff = 1E-4 #correlation bellow treated as zero
        self.expanded = True

        self.Factorization()

    #############################################################################

    def loadReport(self):

        reportFile = os.path.join(self.folder, "report.dat")
        print("Loading random field from file " + reportFile)

        if not os.path.isfile(reportFile): raise ValueError("Error: report file " + reportFile + " not found" )

        self.name = os.path.basename(os.path.normpath(self.folder))
    
        self.cross_correlation_matrix = np.array([[1]]) 

        f = open(reportFile,"r")
        while True:
            line = f.readline()
            if not line: break

            line = line.split()
            if   line[0] =="dimension:": self.dim = int(line[1])
            elif   line[0] =="range_x:": self.x_range = np.array([float(line[1]), float(line[2])])
            elif line[0] =="range_y:": self.y_range = np.array([float(line[1]), float(line[2])])
            elif line[0] =="range_z:": self.z_range = np.array([float(line[1]), float(line[2])])
            elif line[0] =="correlation_length:": self.corr_l = float(line[1])
            elif line[0] =="correlation_function:": self.corr_f = line[1]
            elif line[0] =="number_of_variables:": self.Nvar = int(line[1])
            elif line[0] =="distributions:": 
                k = 1
                self.dist_types = []
                self.dist_params = []
                for var in range(self.Nvar):
                    self.dist_types.append(line[k])
                    pl = int(line[k+1])
                    params = []
                    for l in range(pl):
                        params.append(float(line[k+2+l]))
                    self.dist_params.append(params)
                    k=k+2+pl
            elif line[0] =="cross_correlation:": 
                self.cross_correlation_matrix = np.ones((self.Nvar,self.Nvar))
                k = 1
                for i in range(self.Nvar):
                    for j in range(i+1,self.Nvar):
                        self.cross_correlation_matrix[i,j] = self.cross_correlation_matrix[j,i]  = float(line[k])
                        k = k+1
            elif line[0] =="save_type:": self.filesavetype = line[1]
            elif line[0] =="sampling_type:": self.sampling_type = line[1]
        f.close()


    #############################################################################

    def saveReport(self):
        if not os.path.isdir(self.folder): os.mkdir(self.folder)

        f = open(os.path.join(self.folder,"report.dat"),"w")
        f.write("field_name:\t%s"%self.name + os.linesep)
        f.write("dimension:\t%d"%self.dim + os.linesep)
        f.write("range_x:\t%e\t%e"%(self.x_range[0],self.x_range[1]) + os.linesep)
        f.write("range_y:\t%e\t%e"%(self.y_range[0],self.y_range[1]) + os.linesep)
        f.write("range_z:\t%e\t%e"%(self.z_range[0],self.z_range[1]) + os.linesep)
        f.write("correlation_length:\t%e"%self.corr_l + os.linesep)
        f.write("correlation_function:\t%s"%self.corr_f + os.linesep)
        f.write("number_of_variables:\t%d"%self.Nvar + os.linesep)
        f.write("distributions:")   
        for i in range(self.Nvar):
            f.write("\t%s\t%d"%(self.dist_types[i],len(self.dist_params[i])))
            if self.dist_types[i]=="file": f.write("\t%s"%self.dist_params[i][0])
            else: 
                for j in self.dist_params[i]: f.write("\t%e"%j)
        f.write(os.linesep)
        if (self.Nvar>1):
            f.write("cross_correlation:")
            for i in range(self.Nvar):
                for j in range(i+1,self.Nvar):
                    f.write("\t%e"%self.cross_correlation_matrix[i,j])
            f.write(os.linesep)         
        f.write("save_type:\t%s"%self.filesavetype + os.linesep)
        f.write("sampling_type:\t%s"%self.sampling_type + os.linesep)
        #add other parameters
        f.close()


    #############################################################################

    def saveNumpyFile(self, name, data, fmt="%e"):
        if self.filesavetype == "binary":
            np.save(name, data)
        else:
            np.savetxt(name+".dat", data,fmt=fmt)

    #############################################################################

    def giveNumpyExt(self):
        if self.filesavetype == "binary": return ".npy"
        else: return ".dat"

    #############################################################################

    def loadNumpyFile(self, name):
        if self.filesavetype == "binary":
            return np.load(name + ".npy")
        else:
            return np.loadtxt(name + ".dat")

    #############################################################################

    def CreateDistributionObjects(self):
        self.dist = [] #objects
        self.dist_stat = [] #mean and standard deviation
        for i in range(self.Nvar):
            d1, d2 = crateDistributionObject(self.dist_types[i],self.dist_params[i])
            self.dist.append( d1 )
            self.dist_stat.append( d2 )


    #############################################################################

    def GetCorrelation(self,X2):
        if self.corr_f == "square_exponential": 
            return np.exp(-(X2/self.corr_l2))
        else:
            raise ValueError("Error: autocorrelation function", self.corr_f, "not implemented")

    #############################################################################

    def GenerateCovarianceMatrixExpanded(self,X,lab):
        # print("Generating covariance matrix %s"%lab)
        nx = len(X);
        CX = np.zeros([nx,nx]);
        for i in range(nx):
            dist2 = np.square(X[i]-X);
            corr = self.GetCorrelation(dist2)
            corr[np.where(corr<self.covmat_cutoff)[0]] = 0.
            CX[i,:] = corr
        # print("[DONE]")
        return CX

    #############################################################################
    """
    def GenerateCovarianceMatrixCompactSPARSE(self, nodes):
        print("Generating covariance matrix between grid and nodes")

        n = nodes.shape[0]
        m = self.grid.shape[0]
        row = np.zeros(0).astype(int)
        column = np.zeros(0).astype(int)
        data = np.zeros(0)        
    
        for i in range(n):
            d2 = np.copy(self.grid)
            d2[:,0] -= nodes[i,0]; d2[:,1] -= nodes[i,1]; d2[:,2] -= nodes[i,2];
            d2 = np.sum(np.square(d2),axis=1)
            corr = self.GetCorrelation(d2)
            ind = np.where(corr<self.covmat_cutoff)[0]
            data = np.hstack((data,corr[ind]))
            column = np.hstack((column, ind))
            row = np.hstack((row,column*0+i))
            if i%100==0: print("progress {:2.1%}".format(float(i)/n),)
                    
        nsparse = len(data);
        rowfull = np.zeros(nsparse*self.Nvar); datafull = np.zeros(nsparse*self.Nvar); columnfull = np.zeros(nsparse*self.Nvar)
        p = 0
        for i in range(self.Nvar):
            for j in range(self.Nvar):
                rowfull[p*nsparse:(p+1)*nsparse] = row+j*n
                columnfull[p*nsparse:(p+1)*nsparse+nsparse] = column+i*m
                datafull[p*nsparse:(p+1)*nsparse] = data*self.cross_correlation_matrix[i,j] 
                p += 1


        D = sparse.coo_matrix((datafull, (rowfull, columnfull)), shape=(m*self.Nvar,n*self.Nvar))
        print("[DONE]")
        return D
    """

    #############################################################################

    def GenerateCovarianceMatrixCompact(self, nodes):
        # print("Generating covariance matrix between grid and nodes")

        n = nodes.shape[0]
        m = self.grid.shape[0]
        D0 = np.zeros((m,n))        
    
        for i in range(n):
            d2 = np.zeros(len(self.grid))
            for i in range(dim):
                d2 += np.square(nodes[:,dim])
            corr = self.GetCorrelation(d2)
            ind = np.where(corr<self.covmat_cutoff)[0]
            D0[i,:] = corr
            # if i%100==0: print("progress {:2.1%}".format(float(i)/n),)
                
        D0 = sparse.coo_matrix(D0)
        row = D0.row
        column = D0.col
        data = D0.data
            
        nsparse = len(data);
        rowfull = np.zeros(nsparse*self.Nvar); datafull = np.zeros(nsparse*self.Nvar); columnfull = np.zeros(nsparse*self.Nvar)
        p = 0
        for i in range(self.Nvar):
            for j in range(self.Nvar):
                rowfull[p*nsparse:(p+1)*nsparse] = row+j*n
                columnfull[p*nsparse:(p+1)*nsparse+nsparse] = column+i*m
                datafull[p*nsparse:(p+1)*nsparse] = data*self.cross_correlation_matrix[i,j] 
                p += 1


        D = sparse.coo_matrix((datafull, (rowfull, columnfull)), shape=(m*self.Nvar,n*self.Nvar))
        # print("[DONE]")
        return D

    #############################################################################
    # based on Li HS, Lü ZZ, Yuan XK. Nataf transformation based point estimate method. Chin Sci Bull 2008;53(17):2586–92

    def NatafTransformation(self,var1, var2, corr0, Ngausspoints, p):

        if abs(corr0)<1E-10:return 0.
        if abs(corr0)>1.-1E-10: return 1.

        if Ngausspoints==2:
            Z = [-1,1]
            P = [0.5,0.5]
        elif Ngausspoints==3: 
            Z = [-1.73205080757,0,1.73205080757]
            P = [0.166666666667,0.66666666666,0.166666666667]
        elif Ngausspoints==4: 
            Z = [-2.33441421834,-0.74196378430,0.74196378430,2.33441421834]
            P = [0.045875854768,0.454124145232,0.454124145232,0.045875854768]
        elif Ngausspoints==5: 
            Z = [-2.85697001387,-1.3552617997,0,1.3552617997,2.85697001387]
            P = [0.011257411328,0.222075922006,0.533333333333,0.222075922006,0.011257411328]
        elif Ngausspoints==6: 
            Z = [-3.32425743359,-1.88917587773,-0.616706590154,0.616706590154,1.88917587773,3.32425743359]
            P = [0.00255578440233,0.088615746029,0.408828469542,0.408828469542,0.088615746029,0.00255578440233]
        elif Ngausspoints==7: 
            Z = [-3.75043971768,-2.36675941078,-1.1544053948,0,1.1544053948,2.36675941078,3.75043971768]
            P = [0.000548268858737,0.0307571239681,0.240123178599,0.457142857143,0.240123178599,0.0307571239681,0.000548268858737]
        else: raise ValueError('Number of Guass points out of range 2-7')
    
    
        RN = corr0
        R = 0.
        while abs(RN-R)/corr0>p:
            R = RN; S = 0.; C = np.array([[1., R],[R, 1.]])
            L = np.linalg.cholesky(C)
            for i in range(Ngausspoints):
                for j in range(Ngausspoints):
                    Y = np.dot(L,np.array([Z[i],Z[j]]))                    
                    X = [self.dist[var1].ppf(norm.cdf(Y[0])), self.dist[var2].ppf(norm.cdf(Y[1]))]
                    S = S+P[i]*P[j]*(X[0]-self.dist_stat[var1][0])*(X[1]-self.dist_stat[var2][0])/(self.dist_stat[var1][1] * self.dist_stat[var2][1])
            RN = corr0+R-S
        corr = RN
        return corr

    #############################################################################

    def BackNatafTransformation(self,NG):
        NG2 = np.zeros(NG.shape)
        for var in range(self.Nvar):
            NG2[var] = self.dist[var].ppf(norm.cdf(NG[var]))   
        return NG2


    #############################################################################

    def NatafTransformationInterpolation(self,var):

        if self.dist_types[var]=='Gaussian': 
            return interp1d([-1,1], [-1,1])
        elif os.path.isfile(os.path.join(self.folder, "NatafTransform_%d"%var + self.giveNumpyExt() )):
            # print("Loading Nataf transformation of variable %d"%var)
            x,NR = self.loadNumpyFile(os.path.join(self.folder, "NatafTransform_%d"%var))
            return interp1d(x, NR)

        # print("Nataf transformation of variable %d"%var)
        # Nataf transformation - approximate
        ni = 0
        N = 50
        x = np.linspace(0,1,N)
        NR = np.zeros(N);
        for i in range(N): NR[i] = self.NatafTransformation(var,var, x[i],6,0.005) 

        if 1:
            self.saveNumpyFile(os.path.join(self.folder,'NatafTransform_%d'%var),np.vstack([x,NR]))

            # plt.rcParams.update({'font.size': 18})
            # plt.rcParams.update({'axes.linewidth': 2})
            # plt.rcParams.update({'font.family' : 'serif'})
            # plt.rcParams.update({'font.serif' : 'Times New Roman'})
            # plt.rcParams.update({'text.usetex':True})
    
            fig = plt.figure()
            ax = fig.add_axes([0.15,0.15,0.8,0.8])
            ax.plot([0,1],[0,1],':',color="grey",lw=2)
            ax.plot(x,NR,'-o',color="k",lw=1.5)
            ax.set_xlim([0,1]); ax.set_ylim([0,1])
            #ax.axis("equal")
            ax.set_xlabel("correlation in original space"); ax.set_ylabel("correlation in Gaussian space")            
            fig.savefig(os.path.join(self.folder,'NatafTransform_%d.png'%var))
            #plt.show()
            plt.close(fig)
            
        # print("[DONE]")
        return interp1d(x, NR)

    #############################################################################

    def NatafTransformationCrossCorMatrix(self):

        # print("Nataf transformation of cross correlation matrix")
        self.cross_correlation_matrix_transformed = np.ones((self.Nvar,self.Nvar))
        for i in range(self.Nvar):
            for j in range(i+1,self.Nvar):
                if self.dist_types[i]=='Gaussian' and self.dist_types[j]=='Gaussian':
                    self.cross_correlation_matrix_transformed[i,j] = self.cross_correlation_matrix_transformed[j,i] = self.cross_correlation_matrix[i,j]
                else:
                    self.cross_correlation_matrix_transformed[i,j] = self.cross_correlation_matrix_transformed[j,i] = self.NatafTransformation(i,j, self.cross_correlation_matrix[i,j], 6, 0.005)
        # print("[DONE]")

    #############################################################################

    def FindAllEigenvalues(self,C):
        if C.shape[0]==1: return np.array([1]), np.array([[1]])

        lam, EV =np.linalg.eig(C)
        #lam, EV = sparse.linalg.eigs(C, k=C.shape[0], which='LR', return_eigenvectors=True)
        #lam = np.real(lam)
        #EV = np.real(EV)

        ind = np.flipud(np.argsort(lam))
        lam = lam[ind]
        EV = EV[:,ind]
        return lam, EV 

    #############################################################################

    """
    def FindLargetsEigenvalues(self,C):
        print("searching for eigenvalues and eigenvectors")
        trace = C.shape[0]
        eignum0 = int(trace/8.)
        lam = [0]
        EV = []
        eignum = 0
        while sum(lam)<trace*self.fraction_trace_considered:
            eignum += eignum0
            lam, EV = sparse.linalg.eigs(C, k=eignum, which='LR', return_eigenvectors=True)
            lam = np.real(lam)
            EV = np.real(EV)

        lam = np.array(lam)
        EV = np.array(EV)
        ind = np.flipud(np.argsort(lam))
        lam = lam[ind]
        EV = EV[:,ind]

        t = 0
        sumL = lam[0];
        while sumL<trace*self.fraction_trace_considered:
            t += 1
            sumL += lam[t]
    
        lam = lam[:t]; EV = EV[:,:t]
        #normalize eigenvectors
        #n = len(lam)
        #for i in range(n): 
        #    ind = np.where(abs(EV[:,i])>1./n)[0][0]
        #    if EV[ind,i]<0: EV[:,i] *= -1 

        print("[DONE] used %d eigenvalues, %.3f%% of trace\n"%(len(lam),np.sum(lam)/trace*100.))
        return lam, EV
    """

    #############################################################################

    def Separate2FullEigenvaluesSpatial(self, lambdaX,lambdaY,lambdaZ):

        # print("Combination of eigenmodes")
        nx = len(lambdaX); ny = len(lambdaY); nz = len(lambdaZ)
        n = nx*ny*nz

        lam = np.zeros(nx*ny*nz)
        ijk = np.zeros([nx*ny*nz,3]).astype("int")
        n = 0    
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    lam[i*ny*nz+j*nz+k] = lambdaX[i]*lambdaY[j]*lambdaZ[k]
                    ijk[n,:] = [i,j,k]
                    n += 1

        ind = np.flipud(np.argsort(lam))
        lam = lam[ind]
        ijk = ijk[ind,:]
        
        t = 0
        D = lam/np.sum(lam)
        sumD = 0.
        while sumD<self.fraction_trace_considered:
            sumD = sumD+D[t]
            t += 1

        p = sum(lam[:t])/sum(lam)
        ijk = ijk[:t,:]

        # print("[DONE] used %d eigenvalues, %.3f%% of trace"%(t,p*100))
        return ijk

    #############################################################################

    def Separate2FullEigenvalues(self): 

        #averaged lambdas
        lambdaX = np.zeros(len(self.eigenValsXYZ[0][0]))
        lambdaY = np.zeros(len(self.eigenValsXYZ[0][1]))
        lambdaZ = np.zeros(len(self.eigenValsXYZ[0][2]))

        for var in range(self.Nvar):
            lambdaX = lambdaX + self.eigenValsXYZ[var][0]
            lambdaY = lambdaY + self.eigenValsXYZ[var][1]
            lambdaZ = lambdaZ + self.eigenValsXYZ[var][2]

        lambdaX = lambdaX/self.Nvar
        lambdaY = lambdaY/self.Nvar
        lambdaZ = lambdaZ/self.Nvar
      
# 
        # print("Combination of eigenmodes")
        nx = len(lambdaX); ny = len(lambdaY); nz = len(lambdaZ); nc = len(self.eigenValsC)
        n = nx*ny*nz*nc

        lam = np.zeros(nx*ny*nz*nc)
        ijkc = np.zeros([nx*ny*nz*nc,4]).astype("int")
        n = 0    
        for c in range(nc):
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        lam[c*nx*ny*nz + i*ny*nz+j*nz+k] = lambdaX[i]*lambdaY[j]*lambdaZ[k]*self.eigenValsC[c];
                        ijkc[n,:] = [i,j,k,c]
                        n += 1

        ind = np.flipud(np.argsort(lam))
        lam = lam[ind]
        ijkc = ijkc[ind,:]
        
        t = 0
        D = lam/np.sum(lam)
        sumD = 0.
        while sumD<self.fraction_trace_considered:
            sumD = sumD+D[t]
            t += 1

        p = sum(lam[:t])/sum(lam)
        lam = lam[:t]
        ijkc = ijkc[:t,:]

        # print("[DONE] used %d eigenvalues, %.3f%% of trace"%(t,p*100))
        return ijkc

    #############################################################################
    
    def saveFactorization(self):

        if self.expanded:
            for var in range(self.Nvar):
                self.saveNumpyFile(os.path.join(self.folder,"eigenvalues_X_%d"%var),self.eigenValsXYZ[var][0])
                self.saveNumpyFile(os.path.join(self.folder,"eigenvalues_Y_%d"%var),self.eigenValsXYZ[var][1])
                self.saveNumpyFile(os.path.join(self.folder,"eigenvalues_Z_%d"%var),self.eigenValsXYZ[var][2])
                self.saveNumpyFile(os.path.join(self.folder,"eigenvectors_X_%d"%var),self.eigenVecsXYZ[var][0])
                self.saveNumpyFile(os.path.join(self.folder,"eigenvectors_Y_%d"%var),self.eigenVecsXYZ[var][1])
                self.saveNumpyFile(os.path.join(self.folder,"eigenvectors_Z_%d"%var),self.eigenVecsXYZ[var][2])
            if( self.Nvar>1 ):
                self.saveNumpyFile(os.path.join(self.folder,"eigenvalues_C"),self.eigenValsC)
                self.saveNumpyFile(os.path.join(self.folder,"eigenvectors_C"),self.eigenVecsC)
            self.saveNumpyFile(os.path.join(self.folder,"ijk"),self.ijk,fmt='%d')
        else: 
            self.saveNumpyFile(os.path.join(self.folder,"eigenvalues"),self.lam)
            self.saveNumpyFile(os.path.join(self.folder,"eigenvectors"),self.EV)
            self.saveNumpyFile(os.path.join(self.folder,"gridnodes"),self.grid)




    #############################################################################

    def Factorization(self):       

        #generate grid nodes
        if os.path.isfile(os.path.join(self.folder,"gridnodesX"+self.giveNumpyExt())) and os.path.isfile(os.path.join(self.folder,"gridnodesY"+self.giveNumpyExt())) and os.path.isfile(os.path.join(self.folder,"gridnodesZ"+self.giveNumpyExt())):
            # print("Loading grid")
            GX = self.loadNumpyFile(os.path.join(self.folder,"gridnodesX"))
            GY = self.loadNumpyFile(os.path.join(self.folder,"gridnodesY"))
            GZ = self.loadNumpyFile(os.path.join(self.folder,"gridnodesZ"))
        else:                        
            #CREATE GRID 
            # print("Generating grid")
            xw = (self.x_range_gap[1]-self.x_range_gap[0]); xn = int(np.ceil(xw/self.grid_spacing)+0.5); x0 = self.x_range_gap[0]+(xw-xn*self.grid_spacing)/2.
            yw = (self.y_range_gap[1]-self.y_range_gap[0]); yn = int(np.ceil(yw/self.grid_spacing)+0.5); y0 = self.y_range_gap[0]+(yw-yn*self.grid_spacing)/2.
            zw = (self.z_range_gap[1]-self.z_range_gap[0]); zn = int(np.ceil(zw/self.grid_spacing)+0.5); z0 = self.z_range_gap[0]+(zw-zn*self.grid_spacing)/2.

            GX = np.arange(x0,x0+(xn+0.1)*self.grid_spacing,self.grid_spacing)
            GY = np.arange(y0,y0+(yn+0.1)*self.grid_spacing,self.grid_spacing)
            GZ = np.arange(z0,z0+(zn+0.1)*self.grid_spacing,self.grid_spacing)

            self.saveNumpyFile(os.path.join(self.folder,"gridnodesX"), GX)
            self.saveNumpyFile(os.path.join(self.folder,"gridnodesY"), GY)
            self.saveNumpyFile(os.path.join(self.folder,"gridnodesZ"), GZ)

        X, Y, Z = np.meshgrid(GX,GY,GZ, indexing="ij");
        X = np.reshape(X,[-1,1]); Y = np.reshape(Y,[-1,1]); Z = np.reshape(Z,[-1,1])
        self.grid = np.hstack([X,Y,Z])

        #find eigenvalues and eigenvectors
        if self.expanded:
            allFilesFound = True

            for var in range(self.Nvar):
                if (not os.path.isfile(os.path.join(self.folder,"eigenvalues_X_%d"%var + self.giveNumpyExt()))) or  (not os.path.isfile(os.path.join(self.folder,"eigenvalues_Y_%d"%var + self.giveNumpyExt()))) or (not os.path.isfile(os.path.join(self.folder,"eigenvalues_Z_%d"%var + self.giveNumpyExt()))) or (not os.path.isfile(os.path.join(self.folder,"eigenvectors_X_%d"%var + self.giveNumpyExt()))) or  (not os.path.isfile(os.path.join(self.folder,"eigenvectors_Y_%d"%var + self.giveNumpyExt()))) or (not os.path.isfile(os.path.join(self.folder,"eigenvectors_Z_%d"%var + self.giveNumpyExt()))):
                     allFilesFound = False
            if (self.Nvar>1) and ((not os.path.isfile(os.path.join(self.folder,"eigenvalues_C"+ self.giveNumpyExt()))) or  (not os.path.isfile(os.path.join(self.folder,"eigenvectors_C"+ self.giveNumpyExt())))):
                 allFilesFound = False
            if (not os.path.isfile(os.path.join(self.folder,"ijk"+ self.giveNumpyExt()))):
                 allFilesFound = False

            if allFilesFound:
                # print("Loading factorization")                

                self.eigenValsXYZ = []
                self.eigenVecsXYZ = []
                for var in range(self.Nvar):
                    self.eigenValsXYZ.append([ self.loadNumpyFile(os.path.join(self.folder,"eigenvalues_X_%d"%var)), self.loadNumpyFile(os.path.join(self.folder,"eigenvalues_Y_%d"%var)), self.loadNumpyFile(os.path.join(self.folder,"eigenvalues_Z_%d"%var))])
                    self.eigenVecsXYZ.append([ self.loadNumpyFile(os.path.join(self.folder,"eigenvectors_X_%d"%var)), self.loadNumpyFile(os.path.join(self.folder,"eigenvectors_Y_%d"%var)), self.loadNumpyFile(os.path.join(self.folder,"eigenvectors_Z_%d"%var))])

                if (self.Nvar>1):
                    self.eigenValsC = self.loadNumpyFile(os.path.join(self.folder,"eigenvalues_C"))
                    self.eigenVecsC = self.loadNumpyFile(os.path.join(self.folder,"eigenvectors_C"))
                else:
                    self.eigenValsC = np.ones(1)
                    self.eigenVecsC = np.ones((1,1))  
                    self.cross_correlation_matrix_transformed = self.cross_correlation_matrix

                self.ijk = (self.loadNumpyFile(os.path.join(self.folder,"ijk"))).astype(int)

                self.nataf_interp = []
                for var in range(self.Nvar):
                    self.nataf_interp.append(self.NatafTransformationInterpolation(var))

            else:
            
                # spectral decomposition of cross correlation matrix
                if (self.Nvar>1):
                    self.NatafTransformationCrossCorMatrix()
                    self.eigenValsC, self.eigenVecsC = self.FindAllEigenvalues(self.cross_correlation_matrix_transformed)
                    if(min(self.eigenValsC)<0): raise ValueError("Error: cross correlation matrix is not positive definite")
                    
                else:
                    self.eigenValsC = np.ones(1)
                    self.eigenVecsC = np.ones((1,1))
                    self.cross_correlation_matrix_transformed = self.cross_correlation_matrix

                # compute covariance matrix
                CXraw = self.GenerateCovarianceMatrixExpanded(GX,"X")
                CYraw = self.GenerateCovarianceMatrixExpanded(GY,"Y")
                CZraw = self.GenerateCovarianceMatrixExpanded(GZ,"Z")

                GaussEigenValsX, GaussEigenVecsX = self.FindAllEigenvalues(CXraw)
                GaussEigenValsY, GaussEigenVecsY = self.FindAllEigenvalues(CYraw)
                GaussEigenValsZ, GaussEigenVecsZ = self.FindAllEigenvalues(CZraw)

                self.ijk = self.Separate2FullEigenvaluesSpatial(GaussEigenValsX,GaussEigenValsY,GaussEigenValsZ)

                # spectral decomposition of covariance matrix
                self.nataf_interp = []
                self.eigenValsXYZ = []
                self.eigenVecsXYZ = []

                for var in range(self.Nvar):
                    finterp = self.NatafTransformationInterpolation(var)
                    CX = finterp(CXraw)
                    CY = finterp(CYraw)
                    CZ = finterp(CZraw)
                    self.nataf_interp.append(finterp)
                    eigenValsX, eigenVecsX = self.FindAllEigenvalues(CX)
                    eigenValsY, eigenVecsY = self.FindAllEigenvalues(CY)
                    eigenValsZ, eigenVecsZ = self.FindAllEigenvalues(CZ)
          
                    self.eigenValsXYZ.append([eigenValsX,eigenValsY,eigenValsZ])
                    self.eigenVecsXYZ.append([eigenVecsX,eigenVecsY,eigenVecsZ])                                                                
                                
                self.saveFactorization()
        else:
            pass
            """
            if os.path.isfile(os.path.join(self.folder,"eigenvalues.dat")) and os.path.isfile(os.path.join(self.folder,"eigenvectors.dat")):
                print("Loading factorization")
                self.lam = np.loadtxt(os.path.join(self.folder,"eigenvalues.dat"))
                self.EV = np.loadtxt(os.path.join(self.folder,"eigenvectors.dat"))
            else:
                C = self.GenerateCovarianceMatrixCompact(self.grid)
                self.lam,self.EV = self.FindLargetsEigenvalues(C)
                self.saveFactorization()
            """

        self.Neig = self.ijk.shape[0]*self.Nvar

    #############################################################################

    def saveRandomVariables(self):
        self.saveNumpyFile(os.path.join(self.folder,"random_variables"),self.X)            
    
    #############################################################################
    
    def GenerateRandVariables(self, Nrealz, seed = "none"):

        if seed == "none": 
            seed = np.random.rand()
        np.random.seed(seed)

        self.Nrealz = Nrealz
        if os.path.isfile(os.path.join(self.folder,"random_variables"+self.giveNumpyExt())):
            # print("Loading random variables")
            self.X = self.loadNumpyFile(os.path.join(self.folder,"random_variables"))
            if( not self.X.shape[1] == self.Nrealz):
                raise ValueError("Error: required %d realizations, file %s contains %d realizations"%(self.Nrealz,"random_variables"+self.giveNumpyExt(), self.X.shape[1]))
        else:
            # print("Generating random variables")
            if self.sampling_type=="MC": 
                self.X = np.random.rand(self.Neig,self.Nrealz)
                self.X = norm.ppf(self.X, loc=0., scale=1.)
            elif self.sampling_type=="LHS": 
                self.X = np.zeros([self.Neig,self.Nrealz])
                row = (np.arange(self.Nrealz)+0.5)/self.Nrealz
                for i in range(self.Neig): self.X[i,:]= np.random.permutation(row)
                self.X = norm.ppf(self.X, loc=0., scale=1.)
            elif self.sampling_type=='LHS_RAND':
                self.X = np.zeros([self.Neig,self.Nrealz])
                row = np.arange(self.Nrealz)
                for i in range(self.Neig): self.X[i,:]= (np.random.permutation(row) + np.random.rand(self.Nrealz))/self.Nrealz
                self.X = norm.ppf(self.X, loc=0., scale=1.)
            elif self.sampling_type=='LHS_MED':
                self.X = np.zeros([self.Neig,self.Nrealz])                                
                int_sep = 100
                row = np.zeros(self.Nrealz)
                for i in range(self.Nrealz):     
                    delta = 1./float(self.Nrealz)
                    vals = np.linspace((i+1./(2*int_sep))*delta, (i+1-1./(2*int_sep))*delta, num=int_sep, endpoint=True)
                    vals = norm.ppf(vals)
                    row[i] = np.sum(vals)/int_sep
                for i in range(self.Neig): 
                    self.X[i,:]= np.random.permutation(row)
            elif self.sampling_type=='Freet':
                fin = open(os.path.join(self.folder,"FreetInput.fre"),"w")
                fin.write("[General]\nCountOfSimulations=%d\n"%self.Nrealz)
                fin.write("TypeOfCorrCoef=0\nTypeOfSampling=%d\n\n"%0)
                fin.write("[Category]\nName=X%d\n\n"%self.Nrealz)
                for i in range(len(self.lam)): fin.write("[Variable]\nName=%04d\nDistribution=NORM\nMean=0\nStd=1\nSkew=0\nKurt=0\n\n"%i)
                fin.write("[Category]\nName=Comparative values\n\n") 
                fin.write("[CorrelationMatrix]\n") 
                fin.close()

                curpath = os.getcwd()
                os.chdir(self.folder)
                os.system("FreetInput.fre")
                os.chdir(curpath)

                with open(os.path.join(self.folder,"FreetInput.fre"), "r") as fout: content = fout.read().splitlines()
                conent = np.array(content)
                ind = [i for i,c in enumerate(content) if "[Values]" in content[i]]
                self.X = np.zeros([self.Neig,self.Nrealz])
                for i in  range(self.Neig): 
                    # print(i, self.Neig)
                    self.X[i,:] = [float(val) for val in content[ind[i]+1:ind[i]+1+self.Nrealz]]   
            elif self.sampling_type=='Sobol':

                #this code is needed to install the randtoolbox package
                """
                # import rpy2's package module
                import rpy2.robjects.packages as rpackages

                # import R's utility package
                utils = rpackages.importr('utils')

                # select a mirror for R packages
                utils.chooseCRANmirror(ind=1) # select the first mirror in the list

                packnames = ('randtoolbox')
                from rpy2.robjects.vectors import StrVector
                utils.install_packages(StrVector(packnames))
                """

                # Import packages
                from rpy2.robjects.packages import importr
                randtoolbox = importr('randtoolbox') # Import Functions

                if (self.Neig>1111):
                    raise ValueError("Sorry, R project can generate Sobol sequence only up to dimension 1111, you are requesting dimension ", self.Neig, " - terminating")
                    exit(1)
                Rmatrix = randtoolbox.sobol(self.Nrealz, dim = self.Neig, scrambling = 3, seed = seed)
                self.X = np.array(tuple(Rmatrix)).astype(float).reshape((self.Neig, self.Nrealz))
                self.X = norm.ppf(self.X, loc=0., scale=1.)
            else: 
                raise ValueError("Sampling type ", self.sampling_type, " not implemented - terminating")
                exit(1)
            self.saveRandomVariables()
            # print("[DONE]")



    #############################################################################

    def savePreparation(self):
        for var in range(self.Nvar):
            for i in range(self.PREP[var].shape[1]):
                self.saveNumpyFile(os.path.join(self.folder,"prepData_%d_%03d"%(var,i)),self.PREP[var][:,i])

    #############################################################################

    def precomputeEOLE(self):
       

        allFilesFound = True
        for var in range(self.Nvar):
            for i in range(self.Nrealz): 
                if not os.path.isfile(os.path.join(self.folder,"prepData_%d_%03d"%(var,i) + self.giveNumpyExt() )): 
                    allFilesFound = False
                    break

        if allFilesFound:
            # print('Loading preparation for EOLE')
            self.PREP = []
            for var in range(self.Nvar):            
                self.PREP.append( np.zeros([self.grid.shape[0],self.Nrealz]) )
                for i in range(self.Nrealz): 
                    self.PREP[var][:,i] = self.loadNumpyFile(os.path.join(self.folder,"prepData_%d_%03d"%(var,i))) 
        else: 
        
            # print('Generating preparation for EOLE:')

            nvx = self.eigenVecsXYZ[0][0].shape[0]; nvy = self.eigenVecsXYZ[0][1].shape[0]; nvz = self.eigenVecsXYZ[0][2].shape[0];

            self.PREP = []
            for var in range(self.Nvar): self.PREP.append(np.zeros([nvx*nvy*nvz,self.Nrealz]))

            for var in range(self.Nvar):
                # print("Variable %d out of %d"%(var+1, self.Nvar))

                ntot = self.Neig
                nspatial = len(self.ijk)
                lam = self.eigenValsXYZ[var][0][self.ijk[:,0]]*self.eigenValsXYZ[var][1][self.ijk[:,1]]*self.eigenValsXYZ[var][2][self.ijk[:,2]]
           
                if self.expanded:
                    #collect eigenshapes
                    EV = np.zeros((nspatial,nvx*nvy*nvz))
                    for i in range(nspatial): 
                        EVX = np.einsum("i,j,k->ijk",self.eigenVecsXYZ[var][0][:,self.ijk[i,0]] , self.eigenVecsXYZ[var][1][:,self.ijk[i,1]],self.eigenVecsXYZ[var][2][:,self.ijk[i,2]])
                        EV[i,:] = EVX.flatten()
                
                    #collect eigevalues                    
                    invlamsqrt = 1./np.sqrt(lam)
                    for varx in range(self.Nvar):
                        X = np.einsum("i,ij->ij",invlamsqrt,self.X[nspatial*var:nspatial*(var+1),:])*self.eigenVecsC[varx,var]*np.sqrt(self.eigenValsC[var])
                        self.PREP[varx] += np.einsum("ij,ik->jk",EV,X)
                else:
                    pass
                    """
                    B = np.zeros([self.X.shape[0],self.Nrealz])
                    invlamsqrt = 1./np.sqrt(self.lam)
                    B = np.dot(np.diag(invlamsqrt),self.X)
                    self.PREP = np.zeros([self.Nvar*self.grid.shape[0],self.Nrealz])
                    for i in range(len(self.lam)): 
                        self.PREP += np.dot(np.tile(self.EV[:,i],[self.Nrealz,1]).T,np.diag(B[i,:])) #memory consuming                        
                    """
           
            self.savePreparation()   
    
            # print("[DONE]")


    #############################################################################

    def GetNodes(self, nodefilename= "none"):    
        
        dmin = self.corr_l/6.
        #[specimen,distribution,corrfunction,Nrealzold,ngrid,neig,samtypeold]...
        #= ReadReport(fieldname);

        if nodefilename == "none":
            # print("Generating random nodes",)
            #RANDOM STRUCTURE
            maxiter = 1000; d2min = dmin**2
            iteration = 0
            origin = np.array([self.x_range[0],self.y_range[0],self.z_range[0]])
            size = np.array([self.x_range[1]-self.x_range[0],self.y_range[1]-self.y_range[0],self.z_range[1]-self.z_range[0]])
            nodes = np.expand_dims(np.random.rand(3)*size+origin,axis=0)
            npoints = 1
            while iteration<maxiter:
                node = np.random.rand(3)*size+origin
                d2 = np.copy(nodes)
                d2[:,0] -= node[0]; d2[:,1] -= node[1]; d2[:,2] -= node[2]
                d2 = np.sum(np.square(d2),axis=1)
            
                if min(d2)<d2min: 
                    iteration += 1
                else: 
                    iteration = 0
                    nodes = np.vstack([nodes, node])
                    npoints += 1

            # print("[DONE] generated %d nodes"%nodes.shape[0])

        elif os.path.isfile(nodefilename):
            # print("Loading random nodes from file '%s'"%nodefilename,)
            if (nodefilename.endswith(".npy")):
                nodes = np.load(nodefilename)
            else:
                nodes = np.loadtxt(nodefilename, usecols=range(3))
            # print("[DONE]")
        else: raise ValueError("Error: File with nodes '%s' not found."%snodefilename)
        self.saveNumpyFile(os.path.join(self.folder,'fieldNodes'),nodes);
        return nodes

    #############################################################################

    def GetNodalGaussFieldExpanded(self):
        # print("Generating random field on the GRID")
        nvx = self.eigenVecsXYZ[0][0].shape[0]; 
        nvy = self.eigenVecsXYZ[0][1].shape[0]; 
        nvz = self.eigenVecsXYZ[0][2].shape[0];

        GF = np.zeros([self.Nvar,nvx*nvy*nvz,self.Nrealz])

        for var in range(self.Nvar):
                # print("Variable %d out of %d"%(var+1, self.Nvar))
                ntot = self.Neig
                nspatial = len(self.ijk)
                lam = self.eigenValsXYZ[var][0][self.ijk[:,0]]*self.eigenValsXYZ[var][1][self.ijk[:,1]]*self.eigenValsXYZ[var][2][self.ijk[:,2]]
           
                if self.expanded:
                    #collect eigenshapes
                    EV = np.zeros((nspatial,nvx*nvy*nvz))
                    for i in range(nspatial): 
                        EVX = np.einsum("i,j,k->ijk",self.eigenVecsXYZ[var][0][:,self.ijk[i,0]] , self.eigenVecsXYZ[var][1][:,self.ijk[i,1]],self.eigenVecsXYZ[var][2][:,self.ijk[i,2]])
                        EV[i,:] = EVX.flatten()
                
                    #collect eigevalues                    
                    for varx in range(self.Nvar):
                        X = np.einsum("i,ij->ij",np.sqrt(lam),self.X[nspatial*var:nspatial*(var+1),:])*self.eigenVecsC[varx,var]*np.sqrt(self.eigenValsC[var])
                        GF[varx] += np.einsum("ij,ik->jk",EV,X) 
        # print("[DONE]")
        return GF       

    #############################################################################

    def GetGaussFieldEOLEExpanded(self,nodes):
        # print("generating field on nodes")
        n = nodes.shape[0]
        ng = self.grid.shape[0]
        GF = np.zeros([self.Nvar,n,self.Nrealz])

        for i in range(n):                
            d2 = np.copy(self.grid)
            d2[:,0] -= nodes[i,0]; d2[:,1] -= nodes[i,1]; d2[:,2] -= nodes[i,2]
            d2 = np.sum(np.square(d2),axis=1)
            d2 = self.GetCorrelation(d2)
            d2[np.where(d2<self.covmat_cutoff)[0]]=0.

            for var in range(self.Nvar):
                GF[var,i,:] = np.dot(self.nataf_interp[var](d2),self.PREP[var])
            # if i%1000==0: print("progress {:2.1%}\r".format(float(i)/n),)

        # print("[DONE]")
        return GF

    #############################################################################

    def SaveFieldNodesTXT(self,GF):
        for var in range(self.Nvar):
            self.saveNumpyFile(os.path.join(self.folder,"fieldValues_%d"%var), GF[var])        

    #############################################################################

    def SaveFieldNodesVTKVoronoi(self, realizations = "all"):

        # print('Generating Voronoi tessellation')
        if realizations == "all": realizations = range(self.Nrealz)

        nodes = self.loadNumpyFile(os.path.join(self.folder,"fieldNodes"))

        gap = 0.00001;
        mirrorfactor = 0.8
        supernodes = np.copy(nodes)
        nnodes = nodes.shape[0]    
    
        tomirror = np.where(nodes[:,0]<self.x_range[0]+mirrorfactor*(self.x_range[1]-self.x_range[0]))[0]
        newnodes0 = nodes[tomirror,:]
        newnodes0[:,0] = 2.*self.x_range[0] - newnodes0[:,0] - gap
        tomirror = np.where(nodes[:,0]>self.x_range[1]-mirrorfactor*(self.x_range[1]-self.x_range[0]))[0]
        newnodes1 = nodes[tomirror,:]
        newnodes1[:,0] = 2.*self.x_range[1] - newnodes1[:,0] + gap        
        supernodes = np.vstack((supernodes,newnodes0,newnodes1))
        
        tomirror = np.where(nodes[:,1]<self.y_range[0]+mirrorfactor*(self.y_range[1]-self.y_range[0]))[0]
        newnodes0 = nodes[tomirror,:]
        newnodes0[:,1] = 2.*self.y_range[0] - newnodes0[:,1] - gap
        tomirror = np.where(nodes[:,1]>self.y_range[1]-mirrorfactor*(self.y_range[1]-self.y_range[0]))[0]
        newnodes1 = nodes[tomirror,:]
        newnodes1[:,1] = 2.*self.y_range[1] - newnodes1[:,1] + gap        
        supernodes = np.vstack((supernodes,newnodes0,newnodes1))

        tomirror = np.where(nodes[:,2]<self.z_range[0]+mirrorfactor*(self.z_range[1]-self.z_range[0]))[0]
        newnodes0 = nodes[tomirror,:]
        newnodes0[:,2] = 2.*self.z_range[0] - newnodes0[:,2] - gap
        tomirror = np.where(nodes[:,2]>self.z_range[1]-mirrorfactor*(self.z_range[1]-self.z_range[0]))[0]
        newnodes1 = nodes[tomirror,:]
        newnodes1[:,2] = 2.*self.z_range[1] - newnodes1[:,2] + gap        
        supernodes = np.vstack((supernodes,newnodes0,newnodes1))

        vor = Voronoi(supernodes)

        vnodes = np.zeros(0).astype(int)
        nv = vor.vertices.shape[0]
        indexes = np.zeros(0).astype(int)
        tetras = np.zeros((0,4)).astype(int)

        for i in range(nodes.shape[0]): # range(10):
            cell = np.array(vor.regions[vor.point_region[i]]).astype(int)
            cellnodes = np.vstack((nodes[i,:],vor.vertices[cell,:]))
            tri = Delaunay(cellnodes)
            for k in range(tri.simplices.shape[0]):    
                tetras = np.vstack((tetras,tri.simplices[k,:]+len(vnodes)))
                indexes = np.hstack((indexes,i))
            vnodes = np.hstack((vnodes,np.hstack((nv+i,cell))))
        vnodes, retindices = np.unique(vnodes, return_inverse=True)
        tetras = retindices[tetras]

        # print('Saving VTK Voronoi files')
        for var in range(self.Nvar):
            GF = self.loadNumpyFile(os.path.join(self.folder,"fieldValues_%d"%var))     

            output= open(os.path.join(self.folder,"VTK_Voronoi_%d.vtk"%var),'w')
            output.write("# vtk DataFile Version 2.0\n")
            output.write("Field %s\n"%self.name)
            output.write("ASCII\n\n")

            allvertices = np.vstack((vor.vertices,nodes))    
            
            output.write("DATASET UNSTRUCTURED_GRID\n")
            output.write("POINTS %d double \n"%len(vnodes))
            for i in range(len(vnodes)): output.write("%.5e\t%.5e\t%.5e\n"%(allvertices[vnodes[i],0],allvertices[vnodes[i],1],allvertices[vnodes[i],2]))
            
            output.write("\nCELLS %d %d\n"%(tetras.shape[0],tetras.shape[0]*5))
            for i in range(tetras.shape[0]): output.write("%d\t%d\t%d\t%d\t%d\n"%(4, tetras[i,0], tetras[i,1], tetras[i,2], tetras[i,3]))

            output.write("\nCELL_TYPES %d\n"%tetras.shape[0])
            for i in range(tetras.shape[0]): output.write("10\n")

            output.write("CELL_DATA %d\n"%tetras.shape[0]) 
            for i in realizations:
                field = GF[:,i]
                field = field[indexes]
                output.write("\nSCALARS field_%04d_realz_%04d double 1\n"%(c,i))
                output.write("LOOKUP_TABLE default\n")
                for h in field: output.write("%.5e\n"%h)
            output.close()

    #############################################################################

    def SaveFieldNodesVTKDots(self, realizations = "all"):

        # print('Saving VTK dots files')
        if realizations == "all": realizations = range(self.Nrealz)

        for var in range(self.Nvar):
            GF = self.loadNumpyFile(os.path.join(self.folder,"fieldValues_%d"%var))     
            nodes = self.loadNumpyFile(os.path.join(self.folder,"fieldNodes"))     

            n = nodes.shape[0]

            output= open(os.path.join(self.folder,"VTK_dots_%d.vtk"%var),'w')
            output.write("# vtk DataFile Version 2.0\n")
            output.write("Field %s\n"%self.name)
            output.write("ASCII\n\n")

            output.write("DATASET UNSTRUCTURED_GRID\n")
            output.write("POINTS %d double \n"%n)
            for i in range(n): output.write("%.5e\t%.5e\t%.5e\n"%(nodes[i,0],nodes[i,1],nodes[i,2]))
            
            output.write("\nCELLS %d %d\n"%(n,n*2))
            for i in range(n): output.write("%d\t%d\n"%(1, i))

            output.write("\nCELL_TYPES %d\n"%n)
            for i in range(n): output.write("1\n")

            output.write("CELL_DATA %d\n"%n) 
            for i in realizations:
                field = GF[:,i]
                output.write("\nSCALARS realz_%04d double 1\n"%(i))
                output.write("LOOKUP_TABLE default\n")
                for h in field: output.write("%.5e\n"%h)
            output.close()

    #############################################################################
    
    def SaveGridVTKDots(self, realizations = "all"):
        GF = self.GetNodalGaussFieldExpanded();
        GF = self.BackNatafTransformation(GF); 
        # print('Saving VTK dots files of GRID')
        if realizations == "all": realizations = range(self.Nrealz)

        n = self.grid.shape[0]
        for var in range(self.Nvar):
            output= open(os.path.join(self.folder,"VTK_dots_GRID_%d.vtk"%var),'w')
            output.write("# vtk DataFile Version 2.0\n")
            output.write("Field %s\n"%self.name)
            output.write("ASCII\n\n")

            output.write("DATASET UNSTRUCTURED_GRID\n")
            output.write("POINTS %d double \n"%n)
            for i in range(n): output.write("%.5e\t%.5e\t%.5e\n"%(self.grid[i,0],self.grid[i,1],self.grid[i,2]))
            
            output.write("\nCELLS %d %d\n"%(n,n*2))
            for i in range(n): output.write("%d\t%d\n"%(1, i))

            output.write("\nCELL_TYPES %d\n"%n)
            for i in range(n): output.write("1\n")

            output.write("CELL_DATA %d\n"%n) 
            for i in realizations:
                field = GF[var,:,i]
                output.write("\nSCALARS realz_%04d double 1\n"%(i))
                output.write("LOOKUP_TABLE default\n")
                for h in field: output.write("%.5e\n"%h)
            output.close()

    #############################################################################
    
    def SaveGridVTKVoronoi(self, realizations = "all"):
        GF = self.GetNodalGaussFieldExpanded();
        GF = self.BackNatafTransformation(GF); 

        # print('Saving VTK Voronoi files of GRID')
        if realizations == "all": realizations = range(self.Nrealz)
        
        GX = self.loadNumpyFile(os.path.join(self.folder,"gridnodesX"))
        GY = self.loadNumpyFile(os.path.join(self.folder,"gridnodesY"))
        GZ = self.loadNumpyFile(os.path.join(self.folder,"gridnodesZ"))
        nx = len(GX); ny = len(GY); nz = len(GZ);        

        sx = GX[1]-GX[0]
        GX = np.arange(GX[0]-sx/2.,GX[-1]+sx/1.9,sx)
        sy = GY[1]-GY[0]
        GY = np.arange(GY[0]-sy/2.,GY[-1]+sx/1.9,sy)
        sz = GZ[1]-GZ[0]
        GZ = np.arange(GZ[0]-sz/2.,GZ[-1]+sz/1.9,sz)

        X, Y, Z = np.meshgrid(GX,GY,GZ, indexing="ij");
        X = np.reshape(X,[-1,1]); Y = np.reshape(Y,[-1,1]); Z = np.reshape(Z,[-1,1])
        grid = np.hstack([X,Y,Z])

        n = self.grid.shape[0]
        for var in range(self.Nvar):
            output= open(os.path.join(self.folder,"VTK_Voronoi_GRID_%d.vtk"%var),'w')
            output.write("# vtk DataFile Version 2.0\n")
            output.write("Field %s\n"%self.name)
            output.write("ASCII\n\n")

            output.write("DATASET UNSTRUCTURED_GRID\n")
            output.write("POINTS %d double \n"%(len(grid)))
            for i in range(len(grid)): output.write("%.5e\t%.5e\t%.5e\n"%(grid[i,0],grid[i,1],grid[i,2]))
            
            output.write("\nCELLS %d %d\n"%(n,n*9))
            for i in range(n): 
                ix = int(i/(ny*nz))
                iy = int((i-ix*ny*nz)/nz)
                iz = int(i-ix*ny*nz-nz*iy)
                a = ix*(ny+1)*(nz+1) + iy*(nz+1) + iz

                output.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"%(8, a, a + (ny+1)*(nz+1), a+nz+1 + (ny+1)*(nz+1), a+nz+1, a + 1, a + (ny+1)*(nz+1) + 1, a+nz+1 + (ny+1)*(nz+1) + 1, a+nz+2))

            output.write("\nCELL_TYPES %d\n"%n)
            for i in range(n): output.write("12\n")

            output.write("CELL_DATA %d\n"%n) 
            for i in realizations:
                field = GF[var,:,i]
                output.write("\nSCALARS realz_%04d double 1\n"%(i))
                output.write("LOOKUP_TABLE default\n")
                for h in field: output.write("%.5e\n"%h)
            output.close()
    
    #############################################################################

    def GenerateField(self,nodefile):

        nodes = nodefile
        limits_max = np.max(nodes,axis =0)
        limits_min = np.max(nodes,axis =0)
        nodes_in_box = True
        if (limits_min[0]<self.x_range[0] or limits_max[0]>self.x_range[1]):
            # print("RF x limits: ", self.x_range, "min and max x coordintes of points:", [limits_min[0], limits_max[0]])
            nodes_in_box = False
        if (limits_min[1]<self.y_range[0] or limits_max[1]>self.y_range[1]):
            # print("RF y limits: ", self.y_range, "min and max y coordintes of points:", [limits_min[1], limits_max[1]])
            nodes_in_box = False
        if (limits_min[2]<self.z_range[0] or limits_max[2]>self.z_range[1]):
            # print("RF z limits: ", self.z_range, "min and max z coordintes of points:", [limits_min[2], limits_max[2]])
            nodes_in_box = False
        if not nodes_in_box:
            print("Error in random field: rubmitted points exceed the box in which random field was generated. Termination is strongly suggested as the results will be most likely biased.")
        
 
        GF = self.GetGaussFieldEOLEExpanded(nodes)
        GF = self.BackNatafTransformation(GF)
        return GF


    #############################################################################

    def ErrorEvaluation(self, max_node_num = 5e3):

        print("error evaluation")
            
        GF = []
        nodes = self.loadNumpyFile(os.path.join(self.folder,"fieldNodes"))     

        n = nodes.shape[0]
        param = int(max_node_num) 
        randint = []
        
        print(param,n)
        if n>param: randint = np.random.randint(0,high=n,size=(param))
        else: randint = np.arange(0,n)
        param = len(randint)
        print(param)
        nodesX = nodes[randint,:]

        for var in range(self.Nvar):
            # print("Variable %d out of %d"%(var+1, self.Nvar))            
           
        
            GF.append (self.loadNumpyFile(os.path.join(self.folder,"fieldValues_%d"%var)) )
            
            GFX = GF[var][randint,:]

            # plt.rcParams.update({'font.size': 18})
            # plt.rcParams.update({'axes.linewidth': 2})
            # plt.rcParams.update({'font.family' : 'serif'})
            # plt.rcParams.update({'font.serif' : 'Times New Roman'})
            # plt.rcParams.update({'text.usetex':True})
    
            ################# AUTOCORRELATION
            fig = plt.figure()
            ax = fig.add_axes([0.15,0.15,0.8,0.8])

            #corr,pval = pearsonr(GFX,axis=1)
            corr,pval = spearmanr(GFX,axis=1) 
            for i in range(param):
                d2 = np.copy(nodesX)
                d2[:,0] -= nodesX[i,0]; d2[:,1] -= nodesX[i,1]; d2[:,2] -= nodesX[i,2]
                d2 = np.sum(np.square(d2),axis=1)
                ax.plot(np.sqrt(d2),corr[i,:],'ok',ms=1)
            ax.set_ylabel("Spearman's correlation"); ax.set_xlabel("distance")        
            
            xlim = ax.get_xlim()
            ax.set_xlim([0,xlim[1]])
            x = np.linspace(0,xlim[1],num =100)
            y = self.GetCorrelation(np.square(x))
            ax.plot(x,y,'-r',lw=2)

            fig.savefig(os.path.join(self.folder,'corr_error_%d.png'%var))
            #plt.show()
            plt.close(fig)

            ################# DISTIBUTION
        
            # print(GFX.shape)
            GFXX = np.reshape(GFX,[-1,1])

            hist, bins = np.histogram(GFXX, density=True, bins = int(np.sqrt(len(GFXX)/10.)))    

            fig = plt.figure()
            ax = fig.add_axes([0.15,0.15,0.8,0.8])
            ax.bar((bins[1:]+bins[:-1])/2., hist, color="k", alpha=0.5, width = (bins[1]-bins[0])*0.8)
            xlim = ax.get_xlim()
            x = np.linspace(xlim[0],xlim[1],num =100)
            y = self.dist[var].pdf(x)
            ax.plot(x,y,'-r',lw=2)
            fig.savefig(os.path.join(self.folder,'dist_error_%d.png'%var))
            #plt.show()
            ax.set_ylabel("pdf"); ax.set_xlabel("variable %d"%var)
            plt.close(fig)        

        
        if self.Nvar>1:
            crosscorr = np.zeros((int((self.Nvar**2-self.Nvar)/2+0.5),4))
            t = 0
            for var1 in range(self.Nvar):
                for var2 in range(var1+1,self.Nvar):
                    #correlation = pearsonr(GF[var1].flatten(),GF[var2].flatten())[0]
                    correlation = spearmanr(GF[var1].flatten(),GF[var2].flatten())[0]
                    crosscorr[t,:] = np.array([var1,var2,correlation,self.cross_correlation_matrix[var1,var2]])
                    t += 1
                    
            np.savetxt(os.path.join(self.folder,"cross_correlation_error.dat"),crosscorr)
        
        # print("[DONE]")

# if __name__ == '__main__':

    
#     # RF = RandomField(dist_types=["Uniform"], dist_params=[[1.,4.]], name = "FIELD_Uni", corr_l = 5,filesavetype="text", x_range=[0.,10.],y_range=[0.,10.],z_range=[0.,0.5], sampling_type="MC")
#     # #RF = RandomField(dist_types=["Gaussian"], dist_params=[[3.,0.8]], name = "FIELD_Gaussian", corr_l = 0.5,filesavetype="text")
#     # #RF = RandomField(dist_types=["Gamma"], dist_params=[[1.,0.8]], name = "FIELD_Gamma",corr_l = 0.3, x_range=[0.,1],y_range=[0.,2],z_range=[0.,0.5])
#     # #RF = RandomField(dist_types=["Weibull"], dist_params=[[1,5]], name = "FIELD_Weibull", corr_l = 0.4)

#     # #RF = RandomField(dist_types=["Gaussian","Gaussian"], dist_params=[[3.,0.8],[3.,0.8]], name = "FIELD_Gaussian2x",CC=[[1.,0.8],[0.8,1.0]], corr_l = 0.5,filesavetype="text")
#     # #RF = RandomField(dist_types=["Weibull","Lognormal","Grafted","Gamma"], dist_params=[[1,5],[1,0.25], [1.,0.2, 30,1e-4],[1.,0.8]], name = "FIELD_four_variables",CC=[[1.,0.9, -0.8, 0.1],[0.9,1.0,-0.6,0.],[-0.8,-0.6,1.0,0.2],[0.1,0,0.2,1.0]], corr_l = 0.3, x_range=[0.,1], y_range=[0.,2], z_range=[0.,0.5], sampling_type="LHS_MED")
#     # #RF = RandomField(dist_types=["Grafted","Gamma"], dist_params=[[1.,0.2, 30,1e-4],[1.,0.8]], name = "FIELD_GratedGamma",CC=[[1.,-0.8],[-0.8,1.0]], corr_l = 0.4,filesavetype="npy")
#     # #RF = RandomField(dist_types=["Weibull","Lognormal","Gamma"], dist_params=[[1,5],[1,0.25],[1.,0.8]], name = "FIELD_four_variables",CC=[[1.,0.9, 0.9],[0.9,1.0,0.9],[0.9,0.9,1.0]], corr_l = 0.3, x_range=[0.,1.], y_range=[0.,1.], z_range=[0.,0.1], sampling_type="LHS_MED")
#     # #RF = RandomField(readFromFolder = "FIELD_four_variables")
 

#     # RF.GenerateRandVariables(500, seed = 8)        #generation of random numbers, critical step
#     # #RF.SaveGridVTKDots()                #visualization of GRID field as dots, not needed
#     # #RF.SaveGridVTKVoronoi()             #visualization of GRID field as boxes, not needed
#     # RF.precomputeEOLE()                  #calculation of preparation files, can be used for any LDPM geometry
#     # RF.GenerateField(nodefile="none")    #EOLE projection
#     # #RF.SaveFieldNodesVTKDots()          #visualization as dots, not needed
#     # #RF.SaveFieldNodesVTKVoronoi(realizations = [0,1,2])       #visualization as voronoi cells, not needed
#     # RF.ErrorEvaluation(max_node_num = 5e3)                  #check of correlation, cross correlation and distribution, not needed


#     RF = RandomField(dimension=1, dist_types=["TruncatedGaussian"], dist_params=[[10.,5.,0]], name= "FIELD_1D_TruncNormal", corr_l = 5, filesavetype="text", x_range=[0.,100.],sampling_type="MC")
    
#     RF.GenerateRandVariables(500, seed = 8)
    
#     RF.precomputeEOLE()
    
#     RF.GenerateField(nodefile="none")
    
#     RF.ErrorEvaluation()
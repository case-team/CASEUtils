import numpy as np
import math
import ctypes
import ROOT

(pT_i, eta_i, phi_i) = (0, 1, 2)




def ang_dist(phi1, phi2):
    dphi = phi1 - phi2
    if(dphi < -math.pi):
        dphi += 2.* math.pi
    if(dphi > math.pi):
        dphi -= 2.*math.pi
    return dphi

def raw_moment(jet, eta_order, phi_order):
    return (jet[:,pT_i] * (jet[:,eta_i] ** eta_order) * (jet[:,phi_i] ** phi_order)).sum()
    

def convert_to_pt_eta_phi(jet, jet_conts):
    jet_conv = np.zeros(jet_conts.shape)
    jet_eta = jet[eta_i]
    jet_phi = jet[phi_i]
    pt_sum = 0.
    #Loop thru jet conts and convert from px py pz e to pt eta phi m (is there a way to vectorize this?)
    for i in range(jet_conts.shape[0]):
        if(jet_conts[i,3] <= 0.01): continue #skip 0 energy jets, save time
        vec = ROOT.Math.PxPyPzEVector(jet_conts[i,0], jet_conts[i,1], jet_conts[i,2], jet_conts[i,3])
        jet_conv[i] = [vec.Pt(), vec.Eta() - jet_eta, ang_dist(vec.Phi(), jet_phi), 0.]
        pt_sum += vec.Pt()
    return jet_conv





def pixelate(jet, npix=40, img_width=1.0, rotate = True, norm=True):
    #Augmented version of EnergyFlow package version to add rotation and flipping of jet in pre-processing
    #rotation code based on https://alyssaq.github.io/2015/computing-the-axes-or-orientation-of-a-blob/


    # the image is (img_width x img_width) in size
    pix_width = img_width / npix
    jet_image = np.zeros((npix, npix), dtype = np.float16)

    # remove particles with zero pt
    jet = jet[jet[:,pT_i] > 0]

    # get pt centroid values
    #eta_avg = np.average(jet[:,eta_i], weights=jet[:,pT_i])
    #phi_avg = np.average(jet[:,phi_i], weights=jet[:,pT_i])
    pt_sum = np.sum(jet[:,pT_i])
    m10 = raw_moment(jet, 1,0)
    m01 = raw_moment(jet, 0,1)

    #center image
    eta_avg = m10 / pt_sum
    phi_avg = m01 / pt_sum
    jet[:, eta_i] -= eta_avg
    jet[:, phi_i] -= phi_avg

    if(rotate and jet.shape[0] > 2):
        coords = np.vstack([jet[:,eta_i], jet[:,phi_i]])
        cov = np.cov(coords, aweights=jet[:,pT_i])
        evals, evecs = np.linalg.eig(cov)

        e_max = np.argmax(evals)
        eta_v1, phi_v1 = evecs[:, e_max]  # Eigenvector with largest eigenvalue

        theta = np.arctan((eta_v1)/(phi_v1))
        rotation_mat = np.matrix([[np.cos(theta), -np.sin(theta)],
                          [np.sin(theta), np.cos(theta)]])
        transformed_mat = rotation_mat * coords
        eta_transformed, phi_transformed = transformed_mat.A

        transformed_v1 = rotation_mat * np.vstack([eta_v1, phi_v1])
        t_eta_v1, t_phi_v1 = transformed_v1.A

        #flip so max particle is upper right
        argmax_pt = np.argmax(jet[:,pT_i])
        ptmax_eta, ptmax_phi = eta_transformed[argmax_pt], phi_transformed[argmax_pt]
        if(ptmax_eta < 0): eta_transformed *= -1.
        if(ptmax_phi < 0): phi_transformed *= -1.
    else:
        eta_transformed, phi_transformed = jet[:,eta_i], jet[:,phi_i]

    




    mid_pix = np.floor(npix/2)
    # transition to indices
    eta_indices = mid_pix + np.ceil(eta_transformed/pix_width - 0.5)
    phi_indices = mid_pix + np.ceil(phi_transformed/pix_width - 0.5)

    # delete elements outside of range
    mask = np.ones(jet[:,eta_i].shape).astype(bool)
    mask[eta_indices < 0] = False
    mask[phi_indices < 0] = False
    mask[eta_indices >= npix] = False
    mask[phi_indices >= npix] = False

    #print(zip(eta_indices,phi_indices))
    #print(np.mean(mask))
    eta_indices = eta_indices[mask].astype(int)
    phi_indices = phi_indices[mask].astype(int)

    # construct grayscale image
    for pt,eta,phi in zip(jet[:,pT_i][mask], eta_indices, phi_indices): 
        jet_image[phi, eta] += pt

    # L1-normalize the pt channels of the jet image
    if norm:
        normfactor = np.sum(jet_image)
        if normfactor <= 0.:
            print(jet)
            print(zip(eta_indices,phi_indices))
            print(jet_image)
            raise FloatingPointError('Image had no particles!')
        else: 
            jet_image /= normfactor

    return jet_image

def make_image(jet, jet_conts, npix, img_width= 1.0, rotate = True, norm = True):
    jet_conts = convert_to_pt_eta_phi(jet, jet_conts)
    return np.squeeze(pixelate(jet_conts, npix = npix, img_width = img_width, rotate = rotate, norm = norm))

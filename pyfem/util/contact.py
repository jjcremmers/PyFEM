class Contact:

  def __init__ ( self , props , globdat ):
  
    self.flag = False
    
    if hasattr( props , 'contact' ):
      if hasattr( props.contact , 'type' ):
        self.type = props.contact.type
        self.flag = True
        
  def checkContact ( self , row , val , col , B , globdat ):
  
    print("YY")

'''





def contact_check(row, val, col, B, props, globdat, i_node): 
   
  import numpy as np
  A = globdat.A
  location = globdat.nodes[i_node] + globdat.state.reshape( len(globdat.nodes) , len(globdat.dofs.dofTypes) )[ globdat.dofs.IDmap[i_node] ]

  x = location[0]
  y = location[1]

  if len(location) == 2:
    location_3D = np.zeros( (3) )
    location_3D[:2] = location
  elif len(location) == 3:
    location_3D = location


  if hasattr(props.contact, 'type'):
    contact_type = props.contact.type  #possible types are "inf_line" // "fin_line" // "disc"
  else:
    raise RuntimeError('Please specify contact type to be "inf_line" // "fin_line" or "disc"')
  
  if contact_type != 'inf_line' and contact_type != 'fin_line' and contact_type != 'disc' :
    raise RuntimeError('Please specify contact type to be "inf_line" // "fin_line" or "disc"   Type detected = ' + contact_type )

  ########################################################################################################
  ##### TYPE = INF_LINE ##################################################################################
  ########################################################################################################
  if contact_type == 'inf_line':
    if hasattr(props.contact, 'pnt'):
      pnt = np.array([]) 
      for jj in range(len(props.contact.pnt)):
        if type(props.contact.pnt[jj]) == str:
          pnt = np.append(pnt, eval(props.contact.pnt[jj]))   
        else:
          pnt = np.append(pnt, props.contact.pnt[jj])

    if hasattr(props.contact, 'n'): 
      n = array(props.contact.n)
      n_unit = n / np.linalg.norm(n)
    
    x_ps = location_3D - pnt # arbitrary point on contactline to node location, WATCH OUT FOR THE SIGN OF THIS VALUE HERE
    d = np.dot( x_ps, n_unit )   # distance of potential overshoot
    
    if d <= 0 :
      row, val, col, B = contact_true(row, val, col, B, props, globdat, n_unit, d, i_node)

  ########################################################################################################
  ##### TYPE = FIN_LINE ##################################################################################
  ########################################################################################################
  if contact_type == 'fin_line':
    #if hasattr(props.contact, 'lim'): #format for lim:  lim=[lim_x1, lim_y1, lim_x2, lim_y2]; # aka point 1 and point 2 in between which the line is located
    #  line_lim = array(props.contact.lim)

    if hasattr(props.contact, 'lim'):
      line_lims = np.array([]) 
      for k in range(len(props.contact.lim)):
        if type(props.contact.lim[k]) == str:
          line_lims = np.append(line_lims, eval(props.contact.lim[k]))   
        else:
          line_lims = np.append(line_lims, props.contact.lim[k])

    line_list = line_lims.reshape( int(len(props.contact.lim)/6) , 6 ) # / 6 since each line holds two end points of 3dofs each
    
    for j in range( len(line_list) ):
      line_lim = line_list[j,:]

      #dydx = ( line_lim[1] - line_lim[4] ) / ( line_lim[0] - line_lim[3] )
      #n1 = -1/dydx
      #n = [1.0, n1]
      dir_1 = np.array([ line_lim[3] - line_lim[0], line_lim[4] - line_lim[1], 0])
      n = np.cross( dir_1 , [0,0,1] ) #since second line is always perpendicular due to being 2D, and this line is 3D
      n_unit = n[:2] / np.linalg.norm(n[:2])

      if hasattr(props.contact, 'checkpoint'): #point specified in .pro file to check wheter or not the normal vector is in the correct direction (might be better option than this way)
        correct_pnt = array(props.contact.checkpoint[j*3 : j*3+3])

      x_ps_check = correct_pnt[0:2] - line_lim[0:2] 
      d_check = np.dot( x_ps_check, n_unit )
      if d_check < 0 :
        n_unit *= -1

      x_ps = location_3D[:2] - array([line_lim[0], line_lim[1]]) # arbitrary point on contactline to node location, WATCH OUT FOR THE SIGN OF THIS VALUE HERE
      d = np.dot( x_ps, n_unit )   # distance of potential overshoot
  
      ##### check if point Q btween limits#####
      if d <= 0 :
        over = np.array( [ n_unit[0]*d ,  n_unit[1]*d ] )
        loc_prime = location_3D[:2] - over
        if loc_prime[0] < np.min([line_lim[0], line_lim[3]]) or loc_prime[0] > np.max([line_lim[0], line_lim[3]]): #X-check
          d = 1
        if loc_prime[1] < np.min([line_lim[1], line_lim[4]]) or loc_prime[1] > np.max([line_lim[1], line_lim[4]]): #Y-check
          d = 1

      if d <= 0 :
        row, val, col, B = contact_true(row, val, col, B, props, globdat, n_unit, d, i_node)
      ## HERE ALREADY CHECK WHICH CONTACT LINE IS CLOSEST

      ########################################################################################################################
      # SWITCHING BETWEEN FINITE LINES??? BELOW
      
      else:
        if globdat.iterloop == True:

          #if contact has been made """"VARS FROM NONLINEARSOLVERCONTACT.PY""""
          point_P = location_3D[:2] #huidig punt
          point_add_R = globdat.new_locs[ 2*globdat.dofs.IDmap[i_node] : 2*globdat.dofs.IDmap[i_node]+2 ] #np.array([ 5.0 , 2.0 ]) - locations_3D[i,:2] #proposed location - currect location

          point_Q = line_lim[:2]  #contact line limits
          point_add_S = line_lim[3:5] - line_lim[:2] #direction off contact line
          
          if np.cross( point_add_R , point_add_S ) != 0: #else lines are parallel
            t = np.cross( (point_Q-point_P) , point_add_S ) / np.cross( point_add_R , point_add_S ) # https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
            u = np.cross( (point_Q-point_P) , point_add_R ) / np.cross( point_add_R , point_add_S )
            if t <= 1 and t >= 0 and u <= 1 and u >= 0: #else lines do not intersect within limits of the lines
              globdat.inter = np.vstack( [ globdat.inter, t * point_add_R ] )#np.vstack( [ globdat.inter, point_P + t * point_add_R ] )
              
              #globdat.inter = point_P + t * point_add_R
              IDS = array([2*globdat.dofs.IDmap[i_node], 2*globdat.dofs.IDmap[i_node]+2])
              globdat.IDS = np.vstack( [ globdat.IDS, IDS ] )

              
              inter_unit = ( (t * point_add_R) / np.linalg.norm(t * point_add_R) ) *-1
              inter_d = np.linalg.norm(t * point_add_R)
              row, val, col, B = contact_true(row, val, col, B, props, globdat, inter_unit, inter_d, i_node)
      ######################################################################################################
  
  ########################################################################################################
  ##### TYPE = DISC  ##################################################################################
  ########################################################################################################
  if contact_type == 'disc':

    if hasattr(props.contact, 'center'):
      center_loc = np.array([]) 
      for k in range(len(props.contact.center)):
        if type(props.contact.center[k]) == str:
          center_loc = np.append(center_loc, eval(props.contact.center[k]))   
        else:
          center_loc = np.append(center_loc, props.contact.center[k])

    center_list = center_loc.reshape( int(len(props.contact.center)/3) , 3 ) # / 3 since each center location consisits of x y z coords

    if hasattr(props.contact, 'r'):
      r_list = array( props.contact.r )

    for jj in range( len(center_list) ):
      r = r_list[jj]
      center = center_list[jj,:]

      dista = location_3D - center
      d = np.linalg.norm(dista) - r   

      if d <= 0 :
        n_unit = dista / np.linalg.norm(dista)
        row, val, col, B = contact_true(row, val, col, B, props, globdat, n_unit, d, i_node)  
    
  return row, val, col, B

'''

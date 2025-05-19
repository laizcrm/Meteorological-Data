def shr06(): 
    """
    Plot 0-6 km wind shear vector.

    This function calculates the wind shear between the surface (10 m wind)
    and 6 km altitude using perturbation and base geopotential height fields.
    The wind components at 6 km and at the surface are converted to knots,
    and the shear vector and its magnitude are computed.

    The function then plots:
    - A filled contour of the shear magnitude (in knots).
    - Wind barbs representing the shear vector (difference between 6 km and surface wind).

    The plot is labeled and styled according to predefined parameters.
    """
    # create figure
    plt.figure(figsize=(8,8))
    ph = nc.variables['PH'][time] #perturbation geopotential
    phb = nc.variables['PHB'][time] #base state geopotential
    totalgp = phb + ph # total geopotential
    totalgp = funcs.unstagger(totalgp,'Z') #total geopotential unstaggered
    U = funcs.unstagger(nc.variables['U'][time],'U') # U wind component # UNSTAGGERED
    V = funcs.unstagger(nc.variables['V'][time],'V') # V wind component
    u10kts = conv.ms_to_kts(u10[time]) # sfc wind in kts
    v10kts = conv.ms_to_kts(v10[time]) 
    u6 = funcs.interp_generic(6000, (totalgp/9.81), U) # interp to 6km
    v6 = funcs.interp_generic(6000, (totalgp/9.81), V)
    u6kts = conv.ms_to_kts(u6) # convert 6km wind to kts
    v6kts = conv.ms_to_kts(v6)
    #using 10m wind as sfc wind
    ushr = u6kts - u10kts # calc 0-6 shr in kts
    vshr = v6kts - v10kts
    speed = calc.calc_wspeed(ushr, vshr)
    # plot data
    clevs = np.arange(20,145,5)
    cs = m.contourf(x, y, speed, clevs, cmap=cm.get_cmap('gist_ncar'))
    m.barbs(x[::thin,::thin], y[::thin,::thin], ushr[::thin,::thin], vshr[::thin,::thin],length=opt.barbsize) #plot barbs
    title = '0-6km Shear'
    ftitle = 'shr06-' 
    cblabel = 'kts'
    cbticks = True
    makeplot(cs,title,cblabel,clevs,cbticks,ftitle)

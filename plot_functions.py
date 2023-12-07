
def get_dendrogram_custom(Data, cmap="viridis", savefig="", logscale=True, showplot=False):

    """
    This is a modification of the "get_dendrogram" function from DADApy.

    All rights reverved to the DADApy authors.
    .
    .
    .
    
    Generate a visualisation of the topography computed with ADP.
    The only difference is that it makes optional to show the plot.

    This visualisation fundamentally corresponds to a hierarchy of the clusters built
    with Single Linkage taking as similarity measure the density at the
    border between clusters.
    At difference from classical dendrograms, where all the branches have the same height,
    in this case the height of the branches is proportional to the density of the cluster
    centre.
    To convey more information, the distance in the x-axis between
    clusters is proportional to the population (or its logarithm).

    Args:
        Data: A dadapy data object for which ADP has been already run.
        cmap: (optional) The color map for representing the different clusters,
            the default is "viridis".
        savefig: (str, optional) A string with the name of the file in which the dendrogram
            will be saved. The default is empty, so no file is generated.
        logscale: (bool, optional) Makes the distances in the x-axis between clusters proportional
            to the logarithm of the population of the clusters instead of
            proportional to the population itself. In very unbalanced clusterings,
            it makes the dendrogram more human readable. The default is True.

    Returns:

    """
    # Prepare some auxiliary lists
    e1 = []
    e2 = []
    d12 = []
    L = []
    Li1 = []
    Li2 = []
    Ldis = []
    Fmax = max(Data.log_den)
    Rho_bord_m = np.copy(Data.log_den_bord)
    # Obtain populations of the clusters for fine tunning the x-axis
    pop = np.zeros((Data.N_clusters), dtype=int)
    for i in range(Data.N_clusters):
        pop[i] = len(Data.cluster_indices[i])
        if logscale:
            pop[i] = np.log(pop[i])
    xr = np.sum(pop)
    # Obtain distances in list format from topography
    for i in range(Data.N_clusters - 1):
        for j in range(i + 1, Data.N_clusters):
            dis12 = Fmax - Rho_bord_m[i][j]
            e1.append(i)
            e2.append(j)
            d12.append(dis12)

    # Obtain the dendrogram in form of links
    nlinks = 0
    clnew = Data.N_clusters
    for j in range(Data.N_clusters - 1):
        aa = np.argmin(d12)
        nlinks = nlinks + 1
        L.append(clnew + nlinks)
        Li1.append(e1[aa])
        Li2.append(e2[aa])
        Ldis.append(d12[aa])
        # update distance matrix
        t = 0
        fe = Li1[nlinks - 1]
        fs = Li2[nlinks - 1]
        newname = L[nlinks - 1]
        # list of untouched clusters
        unt = []
        for _ in d12:
            if (e1[t] != fe) & (e1[t] != fs):
                unt.append(e1[t])
            if (e2[t] != fe) & (e2[t] != fs):
                unt.append(e2[t])
            t = t + 1
        myset = set(unt)
        unt = list(myset)
        # Build a new distance matrix
        e1new = []
        e2new = []
        d12new = []
        for j in unt:
            t = 0
            dmin = 9.9e99
            for _ in d12:
                if (e1[t] == j) | (e2[t] == j):
                    if (e1[t] == fe) | (e2[t] == fe) | (e1[t] == fs) | (e2[t] == fs):
                        if d12[t] < dmin:
                            dmin = d12[t]
                t = t + 1
            e1new.append(j)
            e2new.append(newname)
            d12new.append(dmin)

        t = 0
        for _ in d12:
            if (unt.count(e1[t])) & (unt.count(e2[t])):
                e1new.append(e1[t])
                e2new.append(e2[t])
                d12new.append(d12[t])
            t = t + 1

        e1 = e1new
        e2 = e2new
        d12 = d12new

    # Get the order in which the elements should be displayed
    sorted_elements = []
    sorted_elements.append(L[nlinks - 1])

    for jj in range(len(L)):
        j = len(L) - jj - 1
        for i in range(len(sorted_elements)):
            if sorted_elements[i] == L[j]:
                sorted_elements[i] = Li2[j]
                sorted_elements.insert(i, Li1[j])

    add = 0.0
    x = []
    y = []
    label = []
    join_distance = []
    for i in range(len(sorted_elements)):
        label.append(sorted_elements[i])
        j = Data.cluster_centers[label[i]]
        y.append(Data.log_den[j])
        x.append(add + 0.5 * pop[sorted_elements[i]])
        add = add + pop[sorted_elements[i]]
        join_distance.append(add)

    xs = x.copy()
    ys = y.copy()
    labels = label.copy()
    zorder = 0
    for jj in range(len(L)):
        c1 = label.index(Li1[jj])
        c2 = label.index(Li2[jj])
        label.append(L[jj])
        if c1 < len(sorted_elements):
            x.append(join_distance[c1])
        else:
            x.append((x[c1] + x[c2]) / 2.0)
        ynew = Fmax - Ldis[jj]
        y.append(ynew)
        x1 = x[c1]
        y1 = y[c1]
        x2 = x[c2]
        y2 = y[c2]
        zorder = zorder + 1
        plt.plot(
            [x1, x1], [y1, ynew], color="k", linestyle="-", linewidth=2, zorder=zorder
        )
        zorder = zorder + 1
        plt.plot(
            [x2, x2], [y2, ynew], color="k", linestyle="-", linewidth=2, zorder=zorder
        )
        zorder = zorder + 1
        plt.plot(
            [x1, x2], [ynew, ynew], color="k", linestyle="-", linewidth=2, zorder=zorder
        )

    zorder = zorder + 1
    cmal = cm.get_cmap(cmap, Data.N_clusters)
    colors = cmal(np.arange(0, cmal.N))
    plt.scatter(xs, ys, c=labels, s=100, zorder=zorder, cmap=cmap)
    for i in range(Data.N_clusters):
        zorder = zorder + 1
        cc = "k"
        r = colors[labels[i]][0]
        g = colors[labels[i]][1]
        b = colors[labels[i]][2]
        luma = (0.2126 * r + 0.7152 * g + 0.0722 * b) * 255
        if luma < 156:
            cc = "w"
        plt.annotate(
            labels[i],
            (xs[i], ys[i]),
            horizontalalignment="center",
            verticalalignment="center",
            zorder=zorder,
            c=cc,
            weight="bold",
        )
    plt.xlim([-0.02 * xr, xr])
    xname = "Population"
    if logscale:
        xname = "ln(Population)"
        plt.xlim([0, xr])
    plt.xlabel(xname)
    plt.ylabel(r"ln($\rho$)")
    if savefig != "":
        plt.savefig(savefig)
    if (showplot):
        plt.show()

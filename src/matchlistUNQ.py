import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
import os


def mkPlot(x1, y1, x2, y2, lab1, lab2, outname='None', save_dir='./'):

    fig, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(x1, y1, label=f'{lab1}', s=20, color='black')
    ax.scatter(x2, y2, label=f'{lab2}', s=5, color='magenta')

    # ax.set_aspect('equal')

    plt.savefig(os.path.join(
        save_dir, f'{outname}.png'), bbox_inches='tight', dpi=150)
    plt.close()

    return None


def matchIn(master, cat, matchtol, xmas, ymas, xcat, ycat, idcat):
    """
    Input:
    master: an ascii table, n by n, with the relevant information in columns
    cat: the array that will be searched for matches, n by n with relevant
    info in columns
    xmas: the column name where the x-coordinates are listed in the master array
    ymas: the column name where the y-coordinates are listed in the master array
    mmas: the column name where the mag is listed in the master array
    xcat: the column name where the x-coordinates are listed in the "cat"
    (matching) array
    ycat: the column name where the y-coordinates are listed in the "cat"
    (matching) array
    mcat: the column name where the mag is listed in the "cat"
    idcat: column in the matching array with the indices of sources as listed
    in the raw file

    Output:
    master: a shortened version of the master array that was input;
    only sources in the master than have a match in the cat are listed
    match_ids: an len(master) by 1 array that has the indices of the best-
    matching sources from the input cat array (ids come from the id column of
    the cat array; if the cat array was shortened from a longer cat array (to
    which this index list will be applied, the indices should have come from
    that longer cat array)). I'm matching based on index position, although in
    the future, to make things easier (possibly?) use real (kind of random,
    not necessarily listed order) IDs and match based on what the column ID
    says
    """

    matchids_in = np.zeros((len(master), 1))
    master['dist'] = np.zeros((len(master), 1))

    nF = True
    row = 0

    while nF:
        # Finding the difference in the x and y positions between the master
        # array, row by row
        # Checking that the difference in x and y (separately) is smaller than
        # the match tolerance
        # The statements in the brackets return indices for which the
        # conditions are true
        # Then, matchrows is a listing of the rows in the cat array where the
        # sources meet the positional requirements.
        matchrows = cat[(abs(master[f'{xmas}'][row] - cat[f'{xcat}'])
                        <= matchtol) & (abs(master[f'{ymas}'][row] - cat[f'{ycat}'])
                        <= matchtol)]

        # If only one source met the tolerance criteria, the index value for
        # that source is put into the matchids array. It will go in the row
        # corresponding to the row where the master source was.
        if (len(matchrows) == 1):
            matchids_in[row][0] = matchrows[f'{idcat}'][0]
            master['dist'][row] = np.sqrt((master[f'{xmas}'][row]
                                           - matchrows[0][f'{xcat}'])**2
                                          + (master[f'{ymas}'][row]
                                             - matchrows[0][f'{ycat}'])**2)
            row += 1

        # If there is more than one source that meets the criteria,
        # I calculate the distance between all of the match sources
        # and the master source. I put these in an array and find the (row)
        # index of the minimum distance. I proceed to put the cat array
        # source index into the matchids array.
        elif (len(matchrows) > 1):
            distDiff = np.zeros((len(matchrows), 1))
            for dd in range(len(matchrows)):
                distDiff[dd] = np.sqrt((master[f'{xmas}'][row]
                                       - matchrows[f'{xcat}'][dd])**2
                                       + (master[f'{ymas}'][row]
                                       - matchrows[f'{ycat}'][dd])**2)
            small = np.argsort(distDiff)
            idx = small[0]  # the idx of the ff-th distance in distDiff
            matchids_in[row][0] = matchrows[f'{idcat}'][idx]
            master['dist'] = distDiff[idx]
            row += 1

        # If there is nothing that meets the criteria, the master source
        # row is removed, as well as the corresponding row in the matchids
        # array.
        else:
            # print('2',type(master))
            # master = np.delete(master,row,0)
            master.remove_row(row)
            matchids_in = np.delete(matchids_in, row, 0)

        # If the row counter is longer than the length of the master,
        # we've reached the end of the distance tabulations. I do a uniqueness
        # check as if there's a repeat in the match_ids, it means multiple
        # master sources matched with the same cat source. (I'm going
        # methodically through the master list, but not removing best-matching
        # sources from the cat array.)
        if (row >= len(master)):
            matchIDout, udx = np.unique(matchids_in, return_index=True)
            # udx is the array of unique indices. I see if this is less than
            # the master length. If so, I use udx to get the relevant, unique
            # sources.
            if len(udx) < len(master):
                uu = 1
                nF = False
                matchout = matchIDout
                # sorted unique matchID values
                # gets rid of duplicates if a source from "cat" matched multiple "master" sources

            elif len(udx) == len(master):
                uu = 0
                nF = False
                matchout = matchids_in
    return master, matchout, uu


def matchlistID(master, cat, matchtol, xmas, ymas, xcat, ycat, idcat, tag=None, save_dir='./', plot=None):

    m_out, mids, uu = matchIn(master, cat, matchtol,
                              xmas, ymas, xcat, ycat, idcat)
    midx = np.asarray(mids, dtype=int)  # index values to CAT
    # m_out is the trimmed master list, which has a column of indices that relate back to the original master
    # mids are the index values of the cat array that correspond to the cat values that match up to something in m_out
    # if uu==0, it should be a one-to-one; if not, m_out would be longer than mids

    mkPlot(m_out[f'{xmas}'], m_out[f'{ymas}'], cat[midx][f'{xcat}'], cat[midx]
           [f'{ycat}'], lab1='Master1', lab2='Cat1', outname='firstRun', save_dir=save_dir)

    # adding an id column to m_out, because this will be cut down, and we'll need the indices that have a match in shortened cat to pull it out
    m_out['id_u'] = np.arange(0, len(m_out), 1, dtype=int)

    # length master sources with match in catalog (one cat source could match two or more master sources)
    print(f'Len m_out, from master {len(m_out)}')
    # length unique cat values that had matches in master
    print(f'Len mids, from cat {len(mids)}')
    print(f'Len cat, original cat {len(cat)}')  # length original cat

    if uu == 1:
        print('Some overlap in matching indices')
        print(f'Original m_out length: {len(m_out):d}')
        # getting only the unique cat values that had a match in m_out (shortened master); should have 'id' column
        m2 = cat[midx]
        c2 = deepcopy(m_out)  # copying m_out to be put in as new "cat"

        print(f'Len m2, from cat {len(m2)}')
        print(f'Len c2, from master {len(c2)}')

        # Running again with master and cat inverted
        # Each values in the "cat" only appears once
        nF = True
        iter = 0
        while nF:
            # c3 = deepcopy(c2)
            if iter % 2 == 0:
                print(f'ITER {iter:d}')
                print(f'Input length m2, from cat {len(m2)}')
                # m2 is from CAT, c2 is from MAS
                c_out, cids, uu = matchIn(
                    m2, c2, matchtol, xcat, ycat, xmas, ymas, 'id_u')
                cidx = np.asarray(cids, dtype=int)
                print(f'Output length c_out, from cat {len(c_out)}')
                print(f'Output length cids, from master {len(cids)}')
                mkPlot(c_out[f'{xcat}'], c_out[f'{ycat}'], c2[f'{xmas}'], c2[f'{ymas}'], lab1=f'Cat{iter}', lab2=f'Mas{iter}',
                       outname=f'run{iter}', save_dir=save_dir)
            # c_out is the cut-down m2 (copy of cat) (if sources in the cut-down master don't match the cut-down cat anymore)
            # cids are the index values of c2 (copy of m_out) that match with sources in c_out

            else:  # iter is odd,uu==1
                print(f'ITER {iter:d}')
                print(f'Input length m2, from master {len(m2)}')
                # m2 is from MAS, c2 is from CAT
                c_out, cids, uu = matchIn(
                    m2, c2, matchtol, xmas, ymas, xcat, ycat, 'id_u')
                cidx = np.asarray(cids, dtype=int)
            # c_out is cut down m2 (from MAS)
            # cids are index values of c2 (apply to c2 to get array with 'id' values?)
                # should be actual 'id' values from cat
                sCatIdx = np.asarray(c2[cidx][idcat], dtype=int)

                print(f'Output length c_out, from master {len(c_out)}')
                print(f'Output length cids, from cat {len(cids)}')
                mkPlot(c_out[f'{xmas}'], c_out[f'{ymas}'], c2[cidx][f'{xcat}'], c2[cidx][f'{ycat}'], lab1=f'Mas{iter}', lab2=f'Cat{iter}',
                       outname=f'run{iter}', save_dir=save_dir)

            if uu == 0 and iter % 2 == 0:  # uu has gotten to be 0 from the above if/else, even loop; can output and end
                print(f'ITER {iter:d}')
                masterOut = c2[cidx]  # c2 from master here
                # should be the original cat idx values that will match the masterOut list
                match_ids = c_out[f'{idcat}']
                mcidx = np.asarray(match_ids, dtype=int)
                mkPlot(masterOut[f'{xmas}'], masterOut[f'{ymas}'], cat[mcidx][f'{xcat}'], cat[mcidx][f'{ycat}'], lab1=f'Mas{iter}', lab2=f'Cat{iter}',
                       outname=f'end{iter}', save_dir=save_dir)
                nF = False
                # sCat = cat[sCatIdx]
                # sCat['id_u'] = np.arange(0,len(sCat),1,dtype=int)
                # t_out, tids, uu = matchIn(c2[cidx],sCat,matchtol,xmas,ymas,xcat,ycat,'id_u')
                # masterOut = t_out
                # soCat = sCat[tids]
                # match_ids = soCat[idcat]
                # nF = False
            elif uu == 1 and iter % 2 == 0:  # if the problem hasn't been solved after the last iteration, which was even,
                # send to the else loop immediately after while (switch master,cat values below)
                print(f'ITER {iter:d}')
                # individual master values, 2 cuts deep, last copy of c2, which would be from MAS
                m2 = deepcopy(c2[cidx])
                # to be used in matching algorithm with whatever is input as c2
                c_out['id_u'] = np.arange(0, len(c_out), 1, dtype=int)
                # copy of cat, 2 cuts deep with a two cat values that match to one master value
                c2 = deepcopy(c_out)
                iter += 1

            elif uu == 1:  # if still problem, but odd iteration, send to even loop, but switch master and cat
                sCat = cat[sCatIdx]  # shortened cat
                m2 = deepcopy(sCat)
                c2 = c_out
                c2['id_u'] = np.arange(0, len(c_out), 1, dtype=int)
                # m2 = c_out # shortened master
                #
                # c2 = deepcopy(sCat)
                # c2['id_u'] = np.arange(0,len(sCat),1,dtype=int)
                iter += 1

            elif uu == 0:  # if uu==0 but iter is odd;
                masterOut = c_out
                match_ids = sCatIdx
                mcidx = np.asarray(match_ids, dtype=int)
                mkPlot(masterOut[f'{xmas}'], masterOut[f'{ymas}'], cat[mcidx][f'{xcat}'], cat[mcidx][f'{ycat}'], lab1=f'Mas{iter}', lab2=f'Cat{iter}',
                       outname=f'end{iter}', save_dir=save_dir)
                nF = False

                print(f'Length Master Out {len(c_out)}')
                print(f'Length Cat IDs Out {len(match_ids)}')
                # sCat = cat[sCatIdx] # actual cat table with matches to master (has 'id' info)
                # sCat['id_u'] = np.arange(0,len(sCat),1,dtype=int)
                # t_out, tids, uu = matchIn(c2[cidx],sCat,matchtol,xmas,ymas,xcat,ycat,'id_u')
                # masterOut = t_out
                # soCat = sCat[tids]
                # match_ids = soCat[idcat]
                # nF = False
                # m2 = c3[cidx] # individual master values, 2 cuts deep, last copy of c2
                # c_out['id_u'] = np.arange(0,len(c_out),1,dtype=int) #to be used in matching algorithm with whatever is input as c2
                # c2 = deepcopy(c_out) # copy of cat, 2 cuts deep with a two cat values that match to one master value
                # iter += 1

            # print(f'Len cidx {len(cidx)}')
            # print(f'Len c_out {len(c_out)}')

            # if uu==0 and iter%2==0: #c2 should be from MAS here
            #     masterOut = c2[cidx]
            #     match_ids= c_out[f'{idcat}'] # should be the original cat idx values that will match the masterOut list
            #     nF = False
            # # match_ids= cat[c_out[f'{idcat}']]

            # else:
            #     cidx = np.asarray(cids,dtype=int) # the index values to last input c2
            #     m2 = c3[cidx] # individual master values, 2 cuts deep, last copy of c2
            #     c_out['id_u'] = np.arange(0,len(c_out),1,dtype=int) #to be used in matching algorithm with whatever is input as c2
            #     c2 = deepcopy(c_out) # copy of cat, 2 cuts deep with a two cat values that match to one master value
            #     iter += 1

    else:
        masterOut = m_out
        match_ids = mids

    return masterOut, match_ids

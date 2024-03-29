#!/usr/bin/env python

import argparse
import json
import logging
import os
import sys
from bisect import bisect_left

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

logger = logging.getLogger("root")
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)


# ====================================================================
# CLASSES
# ====================================================================
class MyFigure(Figure):
    def __init__(self, *args, **kwargs):
        kwargs["dpi"] = 100
        mpl.rcParams["axes.titlesize"] = 8
        mpl.rcParams["axes.labelsize"] = 8
        mpl.rcParams["xtick.labelsize"] = 8
        mpl.rcParams["ytick.labelsize"] = 8
        mpl.rcParams["legend.title_fontsize"] = 6
        mpl.rcParams["legend.fontsize"] = 6
        super().__init__(*args, **kwargs)

    def savefig(self, *args, **kwargs):
        kwargs["transparent"] = True
        kwargs["bbox_inches"] = "tight"
        kwargs["pad_inches"] = 0.05
        super().savefig(*args, **kwargs)


# ====================================================================
class MyJSON(object):

    """load and store multiple MyJSON class instances in a
    dictionary for a series of times"""

    def __init__(self, indir="output", time=None):
        self.indir = indir
        self.time = time
        self.__set_time_list()
        self.__load_json_data()

    def __set_time_list(self):
        self.__set_time_list_from_files_in_directory()
        if self.time is None:
            pass
        # options below subsample the time list
        elif self.time == "select":
            self.__set_time_list_from_select_times()
        else:
            self.__set_time_list_from_user_input()
        logger.debug("np.shape(time_l)= {}".format(np.shape(self.time_l)))

    def __load_json_data(self):
        self.data_d = {}
        for nn, time in enumerate(self.time_l):
            self.data_d[nn] = self.__get_json_data_from_file(time)

    def get_dict(self, keys=None):
        """get dictionary of data"""
        try:
            return recursive_get(self.data_d, keys)
        except KeyError:
            logger.debug("Key not found: {}".format(keys))
            return None

    def get_time(self, units="yrs"):
        """return time as 2-D numpy array for direct use
        with pyplot, and scale according to units"""
        time_a = np.array(self.time_l, dtype=float)
        if units == "yrs":
            scale = 1.0
        elif units == "kyrs":
            scale = 1.0e-3
        elif units == "gyrs":
            scale = 1.0e-6
        time_a *= scale
        return time_a

    def __set_time_list_from_user_input(self):
        # bit ugly with the parsing here.  Could clean up.
        # if user provides list, do nothing
        if isinstance(self.time, list):
            pass
        else:
            # if user provides comma separated string
            if len(self.time.split(",")) > 1:
                self.time = [int(time) for time in self.time.split(",")]
            # if user provides single value
            else:
                self.time = [int(self.time)]

        # find nearest time in directory to user desired time
        outtime_l = []
        for time in self.time:
            neartime = find_closest(self.time_l, time)
            outtime_l.append(neartime)
        # remove duplicates
        outtime_l = list(dict.fromkeys(outtime_l))
        self.time_l = outtime_l

    def __set_time_list_from_files_in_directory(self):
        """return a time list corresponding to the time information of
        the files in indir"""
        # get times that exist in indir
        try:
            file_l = [
                f
                for f in os.listdir(self.indir)
                if os.path.isfile(os.path.join(self.indir, f))
            ]
        except FileNotFoundError:
            print("no JSON files found in directory: {}".format(self.indir))
            sys.exit(0)
        time_l = [fname for fname in file_l]
        time_l = list(filter(lambda a: a.endswith("json"), time_l))
        time_l = [int(time.split(".json")[0]) for time in time_l]
        # ascending order
        time_l = sorted(time_l, key=int)
        self.time_l = time_l

    def __set_time_list_from_select_times(self):
        """return a list of select times, where subsequent times are
        double the previous time.  This captures the generally trend
        of delayed cooling due to continual outgassing"""
        select_l = [0]
        step = self.time_l[1]  # assume first entry is step size
        desired = step
        while desired < self.time_l[-1]:
            nexttime = find_closest(self.time_l, desired)
            select_l.append(nexttime)
            desired *= 2.0
        select_l.append(self.time_l[-1])  # add last output
        self.time_l = select_l  # update time list

    def __get_json_data_from_file(self, time):
        """load and store json data from file"""
        filename = self.indir + "/{}.json".format(time)
        try:
            json_data = open(filename)
        except FileNotFoundError:
            logger.critical("cannot find file: %s", filename)
            logger.critical("please specify times for which data exists")
            sys.exit(0)
        data_d = json.load(json_data)
        json_data.close()
        return data_d

    def get_values(self, keys, **kwargs):
        # data_l is a list of 1-d np.arrays
        # return None if data not found
        key_l = [0] + keys
        dict_d = self.get_dict(key_l)
        if dict_d is None:
            return None
        data_l = []
        for nn, time in enumerate(self.time_l):
            key_l = [nn] + keys
            data_a = self.__get_values_at_time(key_l, **kwargs)
            data_l.append(data_a)
        # below gives (m,n) (2d array), where m (rows) is
        # time index and n (columns) is data size
        data_a = np.array(data_l)  # , dtype='object' )
        # transpose gives (n,m), which can now be used directly by
        # pyplot to plot data by time (since columns are plotted)
        data_a = data_a.T
        # squeeze dimensions of the array by default
        if kwargs.get("squeeze", True):
            data_a = np.squeeze(data_a)
            # but ensure at least 1d to use consistent indexing for
            # time lists with one entry
            data_a = np.atleast_1d(data_a)
        logger.debug("keys= {}".format("/".join(keys)))
        logger.debug("np.shape= {}".format(np.shape(data_a)))
        return data_a  # 1-d or 2-d np array

    def __get_values_at_time(self, keys, **kwargs):
        """get the scaled values for a particular quantity as a 1d np.array"""
        dict_d = self.get_dict(keys)
        scaling = float(dict_d["scaling"])
        # gets a 1-D array
        values_a = np.asarray(dict_d["values"], dtype=float)
        # but slicing can reduce dimensionality of array, so ensure at
        # least 1-D for consistent processing subsequently
        the_slice = kwargs.get("the_slice", np.index_exp[:])
        values_a = np.atleast_1d(values_a[the_slice])
        logger.debug("np.shape(values_a)= {}".format(np.shape(values_a)))
        scaled_values_a = scaling * values_a
        return scaled_values_a

    # not used here, since it needs LaTeX support and we do not
    # want LaTeX as a dependency
    def get_units(self, keys):
        """get the units (SI) of a particular field"""
        # units do not change with time, so can just access units
        # data from first entry
        key_l = [0] + keys
        dict_d = self.get_dict(key_l)
        units = dict_d["units"]
        units = None if units == "None" else units
        return units


# ====================================================================
def find_closest(myList, myNumber):
    """assumes myList is sorted. Returns closest value to myNumber
    If two numbers are equally close, return the smallest number"""

    # https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value

    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before


# ===================================================================
def plot_entropy(ax, myjson_o):
    """plot entropy"""

    set_xaxis_from_kwargs(ax, myjson_o, xaxis="pressure")
    entropy_b = myjson_o.get_values(["data", "S_b"], squeeze=False)

    # shade melting region grey
    try:
        liquidus = myjson_o.get_values(["data", "liquidus_b"], squeeze=False)
        solidus = myjson_o.get_values(["data", "solidus_b"], squeeze=False)
        ax.fill_between(
            myjson_o.xdata,
            liquidus[:, 0],
            solidus[:, 0],
            facecolor="grey",
            alpha=0.35,
            linewidth=0,
        )

        # dotted lines of constant melt fraction
        for xx in range(0, 11, 2):
            yy_b = xx / 10.0 * (liquidus[:, 0] - solidus[:, 0]) + solidus[:, 0]
            if xx == 0:
                # solidus
                ax.plot(myjson_o.xdata, yy_b, "-", linewidth=0.5, color="black")
            elif xx == 3:
                # typically, the approximate location of the rheological transition
                ax.plot(myjson_o.xdata, yy_b, "-", linewidth=1.0, color="white")
            elif xx == 10:
                # liquidus
                ax.plot(myjson_o.xdata, yy_b, "-", linewidth=0.5, color="black")
            else:
                # dashed constant melt fraction lines
                ax.plot(myjson_o.xdata, yy_b, "--", linewidth=1.0, color="white")
    except TypeError:
        # for single phase, melting curves might not be available
        pass

    for nn, time in enumerate(myjson_o.time_l):
        ax.plot(myjson_o.xdata, entropy_b[:, nn], label=time)

    title = "Entropy"
    set_standard_title(ax, title, "J/kg/K")

    yticks = range(1200, 2801, 400)
    ax.set_yticks(yticks)
    ax.set_ylabel("S", rotation=0)
    ax.yaxis.set_label_coords(-0.15, 0.59)


# ===================================================================
def plot_temperature(ax, myjson_o):
    """plot temperature"""

    set_xaxis_from_kwargs(ax, myjson_o, xaxis="pressure")
    temp_b = myjson_o.get_values(["data", "temp_b"], squeeze=False)

    # plot melting region
    try:
        liquidus = myjson_o.get_values(["data", "liquidus_temp_b"], squeeze=False)
        solidus = myjson_o.get_values(["data", "solidus_temp_b"], squeeze=False)
        ax.fill_between(
            myjson_o.xdata,
            liquidus[:, 0],
            solidus[:, 0],
            facecolor="grey",
            alpha=0.35,
            linewidth=0,
        )
    except TypeError:
        # for single phase systems the phase boundaries might not be exported
        pass

    for nn, time in enumerate(myjson_o.time_l):
        ax.plot(myjson_o.xdata, temp_b[:, nn], label=time)

    title = "Temperature"
    set_standard_title(ax, title, "K")

    yticks = range(1000, 5001, 1000)
    yticks_minor = range(1000, 5001, 500)
    ax.set_yticks(yticks)
    ax.set_yticks(yticks_minor, minor=True)
    ax.set_ylabel("T", rotation=0)
    ax.yaxis.set_label_coords(-0.15, 0.59)


# ===================================================================
def plot_melt_fraction(ax, myjson_o):
    """plot melt fraction (mass fraction of melt)"""

    the_slice = np.index_exp[:]
    set_xaxis_from_kwargs(ax, myjson_o, xaxis="pressure_s", the_slice=the_slice)
    phi = myjson_o.get_values(["data", "phi_s"], the_slice=the_slice, squeeze=False)

    for nn, time in enumerate(myjson_o.time_l):
        ax.plot(myjson_o.xdata, phi[:, nn], label=time)

    title = "Melt fraction"
    set_standard_title(ax, title)

    yticks = np.arange(0, 1.1, 0.2)
    ax.set_ylim(0, 1)
    ax.set_yticks(yticks)
    ax.set_ylabel("Phi", rotation=0)
    ax.yaxis.set_label_coords(-0.125, 0.475)


# ===================================================================
def plot_viscosity(ax, myjson_o):
    """plot viscosity"""

    the_slice = np.index_exp[1:-1]
    set_xaxis_from_kwargs(ax, myjson_o, xaxis="pressure", the_slice=the_slice)
    visc_b = myjson_o.get_values(["data", "visc_b"], the_slice=the_slice, squeeze=False)

    # plot lines for all times
    for nn, time in enumerate(myjson_o.time_l):
        ax.plot(myjson_o.xdata, visc_b[:, nn], label=time)

    title = "Viscosity"
    set_standard_title(ax, title, "Pa s")

    yticks = (1e2, 1e6, 1e12, 1e18, 1e21)
    ax.set_yscale("log")
    ax.set_yticks(yticks)
    ax.set_ylabel("Visc", rotation=0)
    ax.yaxis.set_label_coords(-0.125, 0.67)


# ====================================================================
def set_standard_title(ax, title, units=None):
    """standard title for analysis figures"""

    if units:
        title += ", {}".format(units)
    ax.set_title(title)


# ====================================================================
def set_xaxis_from_kwargs(ax, myjson_o=None, **kwargs):
    """set xaxis based on user input"""

    # default is time on xaxis
    xaxis = kwargs.get("xaxis", "time")

    if xaxis == "time":
        if myjson_o is not None:
            myjson_o.xdata = myjson_o.time_l
        # for magma ocean, almost always will want log10 axis for time
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
        ax.xaxis.set_minor_locator(ticker.LogLocator(base=10, subs=[2, 4, 6, 8]))
        ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
        ax.set_xlabel("Time (yr)")

    elif xaxis == "pressure":
        if myjson_o is not None:
            myjson_o.xdata = myjson_o.get_values(
                ["data", "pressure_b"], squeeze=False, **kwargs
            )[:, 0]
            myjson_o.xdata *= 1.0e-9  # to GPa
        xticks = (0, 50, 100, 135)
        ax.set_xticks(xticks)
        ax.set_xlim(0, 138)
        ax.set_xlabel("P (GPa)")

    elif xaxis == "pressure_s":
        if myjson_o is not None:
            myjson_o.xdata = myjson_o.get_values(
                ["data", "pressure_s"], squeeze=False, **kwargs
            )[:, 0]
            myjson_o.xdata *= 1.0e-9  # to GPa
        xticks = (0, 50, 100, 135)
        ax.set_xticks(xticks)
        ax.set_xlim(0, 138)
        ax.set_xlabel("P (GPa)")

    elif xaxis == "temperature":
        if myjson_o is not None:
            myjson_o.xdata = myjson_o.get_values(["atmosphere", "temperature_surface"])
        max_temp = 2700
        min_temp = 1500
        ticks_temp_major = np.arange(min_temp, max_temp + 1, 250)
        ticks_temp_minor = np.arange(min_temp, max_temp + 1, 125)
        xlabel = "Surface temperature (K)"
        ax.set_xlabel(xlabel)
        ax.set_xlim((max_temp, min_temp))
        ax.set_xticks(ticks_temp_minor, minor=True)
        ax.set_xticks(ticks_temp_major, minor=False)

    elif xaxis == "phi":
        if myjson_o is not None:
            myjson_o.xdata = myjson_o.get_values(
                ["rheological_front_phi", "phi_global"]
            )
        xticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        ax.set_xlim(np.min(xticks), np.max(xticks))
        ax.invert_xaxis()
        ax.xaxis.set_major_locator(ticker.FixedLocator(xticks))
        ax.set_xlabel("Melt fraction")

    elif xaxis == "phi_percent":
        if myjson_o is not None:
            myjson_o.xdata = myjson_o.get_values(
                ["rheological_front_phi", "phi_global"]
            )
            myjson_o.xdata *= 100.0  # to melt percent
        xticks = [0, 20, 40, 60, 80, 100]
        ax.set_xlim(np.min(xticks), np.max(xticks))
        ax.invert_xaxis()
        ax.xaxis.set_major_locator(ticker.FixedLocator(xticks))
        ax.set_xlabel("Melt fraction (\%)")


# ====================================================================
def recursive_get(dict_d, keys=None):
    """function to access nested dictionaries"""

    if not keys:
        return dict_d

    if len(keys) == 1:
        return dict_d[keys[0]]

    return recursive_get(dict_d[keys[0]], keys[1:])


# ===================================================================
def figure_interior(indir="output", time=None):
    """plot standard interior fields"""

    myjson_o = MyJSON(indir, time)

    fig = plt.figure(FigureClass=MyFigure, figsize=(4.7747, 4.7747))
    axs = fig.subplots(2, 2)
    fig.subplots_adjust(wspace=0.3, hspace=0.4)

    plot_entropy(axs[0][0], myjson_o)
    plot_temperature(axs[0][1], myjson_o)
    plot_melt_fraction(axs[1][0], myjson_o)
    plot_viscosity(axs[1][1], myjson_o)

    # legend
    handles, labels = axs[0][1].get_legend_handles_labels()
    # get labels from time_l above and format
    labels = ("{:.2e}".format(time) for time in myjson_o.time_l)
    axs[1][1].legend(handles, labels, ncol=2, title="Time (yrs)")

    fig.savefig("interior.pdf")


# ===================================================================
def main():
    # arguments (run with -h to summarize)
    parser = argparse.ArgumentParser(description="SPIDER lite plotting script")
    parser.add_argument(
        "-t",
        "--times",
        default=None,
        type=str,
        help="Comma-separated list of times (default: all times)",
    )
    parser.add_argument(
        "-d",
        "--directory",
        default="output",
        type=str,
        help="Directory of output SPIDER JSON files (default: output)",
    )

    args = parser.parse_args()
    figure_interior(indir=args.directory, time=args.times)

    plt.show()


# ===================================================================
if __name__ == "__main__":
    main()

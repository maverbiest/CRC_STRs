#!/usr/bin/env python3

import re
from collections import OrderedDict

from tral.repeat import repeat
from tral.repeat_list import repeat_list
from tral.sequence import repeat_detection_io, repeat_detection_run


class DetectorPhobosMax(repeat_detection_run.DetectorPhobos):        
    def phobos_get_repeats(self, infile):
        """ Read repeats from a PHOBOS output file stream successively.
        Read repeats from a PHOBOS output file stream successively.
        Postcondition: infile points to EOF.
        Args:
            infile (file stream): File stream from PHOBOS output.
        Returns:
            (Repeat): A generator function is returned that, when called in a loop,
            yields one repeat per iteration.
        .. todo:: Show PHOBOS output syntax.
        """

        pattern_begin = re.compile(r"(\d+) :\s+\d")
        pattern_seq = re.compile(r"(?<![0-9])([\-ACGT]+)")
        pattern_indels = re.compile(r"(?:\(([0-9]+,[DI])\))")

        # Our possible parser states:
        #
        # state 1: Find TR begin in sequence
        # state 2: Find insertions
        # state 3: Find actual TR sequence
        # state 4: Find ideal TR sequence, extract repeat units and yield

        state = 1
        for i, line in enumerate(infile):
            # print("Line %d: %s", i, line[0:-1])
            if 1 == state:  # Find TR offset
                search = pattern_begin.search(line)
                if search and search.groups()[0] is not None:
                    # print(" *(1->2) Found tandem repeat begin")     
                    # Find TR begin in sequence, init RepeatRegion               
                    region = repeat_detection_io.RepeatRegion()
                    region.begin = int(search.groups()[0])
                    region.msa = []
                    # Find TR concensus unit
                    unit_match = pattern_seq.search(line)
                    unit = unit_match.groups()[0]
                    # Update state
                    state = 2

            elif 2 == state:  # Find potential insertions
                # Extract region insertions and deletions as list of tuples ([0-9]+, [DI])
                indels = [(int(i.split(",")[0]), i.split(",")[1]) for i in pattern_indels.findall(line)]
                del_counter = 0
                insertion_sites = []
                # The position of an insertion needs to be incremented once for each deletion that
                ## occurs to the left of it in the repeat region
                for indel in indels:
                    if indel[1] == "D":
                        del_counter += 1
                        continue
                    elif indel[1] == "I":
                        insertion_sites.append(int(indel[0]) - 1 + del_counter)
                state = 3

            elif 3 == state:  # Find actual TR sequence
                match = pattern_seq.search(line)
                if match and match.groups()[0] is not None:
                    # print(" *(2->3) Found actual repeat region")
                    actual_tr = match.groups()[0]
                    state = 4

            elif 4 == state:  # Find ideal TR sequence, extract repeat units and yield
                match = pattern_seq.search(line)
                if match and match.groups()[0] is not None:
                    # print(" *(3->3) Found ideal repeat region")
                    ideal_tr = match.groups()[0]
                    
                    # get locations of unit start and end sites, as well as gap position relative to unit starts
                    unit_locations_raw, gaps = self.map_unit_boundaries(ideal_tr, unit, insertion_sites)
                    # remove out of range units, trim border-overlapping units
                    trimmed_unit_locations = self.trim_unit_locations(unit_locations_raw, len(unit), len(actual_tr))                    
                    # extract units from TR to construct msa      
                    for i in trimmed_unit_locations:
                        trimmed_unit_locations[i]["unit"] = actual_tr[trimmed_unit_locations[i]["begin"]:trimmed_unit_locations[i]["end"]+1]                    
                    
                    # Get first unit, fill out by adding "-" to start of unit until it is as long as the consensus unit
                    lead_unit = list(trimmed_unit_locations.items())[0][1]
                    if len(lead_unit["unit"]) != len(unit):
                        difference = len(unit) - len(lead_unit["unit"]) 
                        lead_unit["unit"] = "-" * difference + lead_unit["unit"]
                    # Get last unit, fill out by adding "-" to end of unit until it is as long as the consensus unit
                    tail_unit = list(trimmed_unit_locations.items())[-1][1]
                    if len(tail_unit["unit"]) != len(unit):
                        difference = len(unit) - len(tail_unit["unit"]) 
                        tail_unit["unit"] = tail_unit["unit"] + "-" * difference 

                    # insert "-" at the gap locations in each unit
                    ## make sure to not add "-" if the insert causing the gap is in the current unit
                    if gaps:
                        for i in trimmed_unit_locations: # will return keys
                            counter = 0
                            # for gap in unit_gaps:
                            for gap in gaps:
                                if gap in trimmed_unit_locations[i]["inserts"]:
                                    # unit contains gap-causing insert, do not add "-" to this unit
                                    counter += 1
                                    continue
                                # correct for potential insertion of earlier gap in unit
                                gap += counter
                                counter += 1
                                # insert "-" at proper position in unit
                                trimmed_unit_locations[i]["unit"] = trimmed_unit_locations[i]["unit"][:gap] + \
                                    "-" + trimmed_unit_locations[i]["unit"][gap:]
                    # generate msa list from OrderedDict, msa_list is list of strings representing the repeat units
                    msa_list = [trimmed_unit_locations[i]["unit"] for i in trimmed_unit_locations]

                    # print(" *(3->1) repeat region finished, yielding.")
                    state = 1                
                    if len(msa_list) >= 2:
                        region.msa = msa_list
                        yield region
                    else:                        
                        print("phobos: Msa too short %s", str(region.msa))

    def map_unit_boundaries(self, ideal_tr, unit, inserts):   
        """ Map rough unit begin and start sites on repeat region
        Start by finding a match for the consensus TR unit in ideal TR representation provided by PHOBOS.
        Starting from this unit, generate a left and a right flank of units until whole repeat region is covered. Merge 
        these three components and return as one OrderedDict
        This will likely contain units that extend beyond the region begin and end, which will be trimmed later
        using trim_unit_locations(). Unit begin and start sites for left and right flank will be corrected for
        insertions using correct_inserts_left() and correct_inserts_right(), respectively.

        Parameters
        ideal_tr (str):     String representation of the ideal repeat region, extracted from PHOBOS output
        unit (str):         String representation of consensus TR unit, extracted from PHOBOS output
        inserts (list(int)):     
                            List of integers specifying where in the repeat region insertions have occurred,
                            extracted from PHOBOS output (preprocessed to account for deletions prior to inserts)

        Returns
        merged_odict (OrderedDict):
                            An ordered dictionary where the keys are (ascending) integers, and values are dictionaries
                            describing the begin, end and insertion sites for each unit
        gaps (list(int)):   A list of gap positions (relative to unit start site) that will be needed to insert gaps
                            in proper places later
        """     
        region_l = len(ideal_tr)
        unit_l = len(unit)
        i = 0
        first_match = ""
        # Slide consensus unit over ideal sequence until match is found. Extract the boundaries of this match
        ## and construct 'first_match' tuple
        while i <= region_l:
            if ideal_tr[i:i + unit_l] == unit:
                first_match = (i, i + unit_l - 1)
                break
            i += 1
        if not first_match:
            raise ValueError("No match for unit '{}' found in ideal TR region".format(unit))                

        # initialize and populate left flank as (begin, end)
        ## keep adding units to the left of first_match until out of range of region
        left_flank = [(first_match[0] - unit_l, first_match[1] - unit_l)]
        while left_flank[0][1] >= 0:
            new_unit =(
                left_flank[0][0] - unit_l,
                left_flank[0][1] - unit_l
            )
            # insert unit tuples at begin of list so they remain in order of appearance in repeat region
            left_flank.insert(0, new_unit)

        # initialize and populate right flank
        ## keep adding units to the right of first_match until out of range of region
        right_flank = [(first_match[0] + unit_l, first_match[1] + unit_l)]
        while right_flank[-1][0] <= region_l - 1:
                    new_unit = (
                        right_flank[-1][0] + unit_l, 
                        right_flank[-1][1] + unit_l
                    )
                    right_flank.append(new_unit)

        # Now that all (begin, end) tuples have been generated, an OrderedDict is made for both flanks and for the first match.
        ## Keys will be ascending integers and values will be dicts {begin: int, end: int, inserts: list}
        key_counter = 1
        left_flank_odict = OrderedDict()
        for i in left_flank:
            left_flank_odict[key_counter] = {
                "begin": i[0],
                "end": i[1],
                "inserts": []
            }
            key_counter += 1

        first_match_odict = OrderedDict([(key_counter, {"begin": first_match[0], "end": first_match[1], "inserts": []})])
        key_counter += 1

        right_flank_odict = OrderedDict()
        for i in right_flank:
            right_flank_odict[key_counter] = {
                "begin": i[0],
                "end": i[1],
                "inserts": []
            }
            key_counter += 1
        
        # Modify unit begin and end sites to correct for insertions
        ## Also extract gap postions relative to unit begin for later use
        all_gaps = None
        if inserts:
            left_flank_odict, left_gaps = self.correct_inserts_left(left_flank_odict, inserts)
            right_flank_odict, right_gaps = self.correct_inserts_right(right_flank_odict, inserts)
            all_gaps = left_gaps.union(right_gaps)
            all_gaps = sorted(list(all_gaps))

        # Merge the three OrderedDicts and return
        ## Note: out of range units are still included here, they will be removed/trimmed later
        merged_odict = OrderedDict(list(left_flank_odict.items()) + list(first_match_odict.items()) + list(right_flank_odict.items()))
        return merged_odict, all_gaps

    def correct_inserts_left(self, left_flank, insertion_sites):  
        """ Correct for insertion sites occurring in the left flank of TR region.
        When an insertion occurs, the boundaries of units in the left flank have to be corrected and shifted
        toward the begin of the repeat region

        Parameters
        left_flank (OrderedDict):      
                                Ordered dictionary with ascending integers as keys and dictionaries describing repeat units as values.
                                e.g. {1:{begin: int, end: int, insertions: list(int)}} Unit begin and end positions are to be
                                corrected for insertions
        insertion_sites (list):
                                List of integers specifying where in the repeat region insertions have occurred

        Returns
        left_flank (OrderedDict):
                                Same Ordered dictionary as input, but with the corrected begin and end position of TR units
        gap_positions (set):    A set of relative gap positions. Relative meaning that gap position are counted from unit begin site,
                                not region begin site
        """      
        gap_positions = set()
        for insert in insertion_sites:
            if insert > list(left_flank.items())[-1][1]["end"]: # insert is in right flank, skip                
                continue
            for i in range(list(left_flank.items())[-1][0], 0, -1):
                if insert >= left_flank[i]["begin"] and insert <= left_flank[i]["end"]:
                    # Insert in unit: push begin back one position, update gaps                    
                    left_flank[i]["begin"] -= 1                       
                    insert_site = insert - left_flank[i]["begin"]                                     
                    left_flank[i]["inserts"].append(insert_site)
                    gap_positions.add(insert_site)                    
                elif insert > left_flank[i]["end"]:
                    # Insert ahead of unit, push both begin and end back 
                    left_flank[i]["begin"] -= 1
                    left_flank[i]["end"] -= 1      
        return left_flank, gap_positions

    def correct_inserts_right(self, right_flank, insertion_sites):
        """ Correct for insertion sites occurring in the right flank of TR region.
        When an insertion occurs, the boundaries of units in the right flank have to be corrected and shifted
        toward the end of the repeat region

        Parameters
        right_flank (OrderedDict):      
                                Ordered dictionary with ascending integers as keys and dictionaries describing repeat units as values.
                                e.g. {1:{begin: int, end: int, insertions: list(int)}} Unit begin and end positions are to be
                                corrected for insertions
        insertion_sites (list):
                                List of integers specifying where in the repeat region insertions have occurred

        Returns
        right_flank (OrderedDict):
                                Same Ordered dictionary as input, but with the corrected begin and end position of TR units
        gap_positions (set):    A set of relative gap positions. Relative meaning that gap position are counted from unit begin site,
                                not region begin site
        """
        gap_positions = set()
        for insert in insertion_sites:
            if insert < list(right_flank.items())[0][1]["begin"]: # insert is in left flank, skip                
                continue
            for i in range(list(right_flank.items())[0][0], list(right_flank.items())[-1][0] + 1):
                if insert >= right_flank[i]["begin"] and insert <= right_flank[i]["end"]:
                    # Insert in unit: push end up one position, update gaps
                    right_flank[i]["end"] += 1                        
                    insert_site = insert - right_flank[i]["begin"]                                    
                    right_flank[i]["inserts"].append(insert_site)
                    gap_positions.add(insert_site)                    
                elif insert < right_flank[i]["begin"]:
                    # Insert prior to unit, push both begin and end up
                    right_flank[i]["begin"] += 1
                    right_flank[i]["end"] += 1             
        return right_flank, gap_positions

    def trim_unit_locations(self, units_odict, unit_length, region_length):
        """ Correct units that overlap repeat region begin or end (e.g. (-3, 4) -> (0, 4)) 
        and trim units that fall out of the repeat region range (e.g. remove (-10, -3))

        Parameters
        units_odict (OrderedDict):       
                                Ordered dictionary with ascending integers as keys and dictionaries describing repeat units as values.
                                e.g. {1:{begin: int, end: int, insertions: list(int)}}
                                Note: unit boundaries should have been corrected for insertions by the
                                appropriate functions previously
        unit_length (int):      Length of the consensus repeat unit
        region_length (int):    Length of the repeat region

        Returns
        units_odict (OrderedDict):
                                Processed input OrderedDict: Out of range units have been removed and 
                                region boundary overlaps have been corrected
        """
        for i in range(1, len(units_odict) + 1):
            if units_odict[i]["end"] < 0:
                del units_odict[i]
            elif units_odict[i]["begin"] >= region_length:
                del units_odict[i]
            elif units_odict[i]["begin"] < 0:
                units_odict[i]["begin"] = 0
            elif units_odict[i]["end"] > region_length:
                units_odict[i]["end"] = region_length
        return units_odict

def main():    
    detector = DetectorPhobosMax()
    detector.config.valopts["--reportUnit"] = 0
    detector.config.valopts["--maxUnitLen"] = 15

    input_file = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/test_phobos/phobos_parsing_test.txt"
    # input_file = "/cfs/earth/scratch/verb/projects/CRC_STRs/results/test/test_phobos/avg_prot_phobos_asis.txt"
    with open(input_file, "r") as f:
        tr_list = list(detector.phobos_get_repeats(infile=f))

    repeats = [repeat.Repeat(i.msa, i.begin) for i in tr_list]       

    repeats_l = repeat_list.RepeatList(repeats)
    for i in repeats_l.repeats:
        print(i)


if __name__ == "__main__":
    main()

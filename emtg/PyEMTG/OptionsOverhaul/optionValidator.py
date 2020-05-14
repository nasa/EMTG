# EMTG: Evolutionary Mission Trajectory Generator
# An open-source global optimization tool for preliminary mission design
# Provided by NASA Goddard Space Flight Center
#
# Copyright (c) 2013 - 2020 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Other Rights Reserved.
#
# Licensed under the NASA Open Source License (the "License"); 
# You may not use this file except in compliance with the License. 
# You may obtain a copy of the License at:
# https://opensource.org/licenses/NASA-1.3
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
# express or implied.   See the License for the specific language
# governing permissions and limitations under the License.

#validator for EMTGv9 journeyoptions/missionoptions entries
#Jacob Englander 1/10/2019

def validate(option):
    if 'name' not in option:
        raise 'missing name field'
    elif option['name'] == '':
        raise 'empty name field'


    if 'dataType' not in option:
        raise 'missing dataType in option "' + option['name'] + '"'
    elif option['dataType'] == '':
        raise 'empty dataType field in option "' + option['name'] + '"'


    if 'defaultValue' not in option:
        raise 'missing defaultValue in option "' + option['name'] + '"'
    elif option['defaultValue'] == '':
        raise 'empty defaultValue field in option "' + option['name'] + '"'

    if option['dataType'] != 'std::string':
        if 'lowerBound' not in option:
            raise 'missing lowerBound in option "' + option['name'] + '"'
        elif option['lowerBound'] == '':
            raise 'empty lowerBound field in option "' + option['name'] + '"'

        if 'upperBound' not in option:
            raise 'missing upperBound in option "' + option['name'] + '"'
        elif option['upperBound'] == '':
            raise 'empty upperBound field in option "' + option['name'] + '"'

    if 'std::vector' in option['dataType']:
        if 'length' not in option:
            raise 'missing length in option "' + option['name'] + '"'
        elif option['length'] == '':
            raise 'empty length field in option "' + option['name'] + '"'


    if 'comment' not in option:
        raise 'missing comment in option "' + option['name'] + '"'
    elif option['comment'] == '':
        raise 'empty comment field in option "' + option['name'] + '"'


    if 'description' not in option:
        raise 'missing description in option "' + option['name'] + '"'
    elif option['description'] == '':
        raise 'empty description field in option "' + option['name'] + '"'

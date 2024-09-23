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

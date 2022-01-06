


def getReads(read_type, SRR_, ext):
    if read_type == "paired_end_read":
        read_1 = "{}_1{}".format(SRR_, ext)
        read_2 = "{}_2{}".format(SRR_, ext)
    elif(read_type == "single_end_read"):
        read_1 = "{}_1{}".format(SRR_, ext)
        read_2 = ""
    else:
        read_1 = ""
        read_2 = ""
    return read_1, read_2
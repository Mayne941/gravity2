import numpy as np


def make_labels(dendro, labels):
    Taxo2ClassDict = {
        TaxoLabel: TaxoGrouping for TaxoLabel, TaxoGrouping in labels}

    for Clade in dendro.find_clades(terminal=True):
        Clade.name = Taxo2ClassDict[Clade.name]

    classes, lines = [], [-1]
    TerminalNodeList = [
        TerminalNode for TerminalNode in dendro.get_terminals()]
    while len(TerminalNodeList) != 0:
        FarLeftNode = TerminalNodeList[0]
        for Clade in ([dendro]+dendro.get_path(FarLeftNode)):
            DescendantNodeList = Clade.get_terminals()
            DescendantClassLabelList = list(
                set([c.name for c in DescendantNodeList]))
            if len(DescendantClassLabelList) == 1:
                classes.append(DescendantClassLabelList[0])
                lines.append(lines[-1]+len(DescendantNodeList))
                TerminalNodeList = TerminalNodeList[len(
                    DescendantNodeList):]
                break

    return np.array(classes), np.array(lines) + 0.5

def split_labels(OrderedTaxoLabelList):
    '''Clean join and split labels between X and Y axes (by species and virus name)'''
    ClassLabelList_minor = []
    for label in OrderedTaxoLabelList:
        label_split = label.split("_")[1:]
        label_final = []
        for word in label_split:
            if "-" in word:
                word = f"---{word.replace('-', ' ')}"
            label_final.append(word)
        ClassLabelList_minor.append(f"{' '.join(label_final)}")

    ClassLabelList_x = [i.split("---")[-1] for i in ClassLabelList_minor]
    ClassLabelList_y = [i.split("---")[0] for i in ClassLabelList_minor]
    return ClassLabelList_x, ClassLabelList_y

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

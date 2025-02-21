import pytest
import matplotlib.pyplot as plt
from app.utils.heatmap_params import get_hmap_params
from app.utils.heatmap_params import construct_hmap_lines

def test_get_hmap_params():
    n_viruses = 50
    n_pphmms = 99
    is_square = True
    virus_dendrogram = None
    do_logscale = False

    hmap_params, fig, ax_dendrogram, ax_Heatmap = get_hmap_params(n_viruses, n_pphmms, is_square, virus_dendrogram, do_logscale)

    assert isinstance(hmap_params, dict)
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax_dendrogram, plt.Axes)
    assert isinstance(ax_Heatmap, plt.Axes)

def test_construct_hmap_lines():
    fig, ax_Heatmap = plt.subplots()
    len_x = 50
    LineList_major = [10, 20, 30, 40, 50, 60]
    LineList_minor = [5, 15, 25, 35, 45, 55]
    hmap_params = {
        'Outer_margin': 0.5,
        'FontSize': 6,
        'dpi': 100,
        'Heatmap_width': 12.0,
        'Heatmap_height': 12.0,
        'TaxoLable_space': 1.00,
        'CBar_Heatmap_gap': 0.05,
        'CBar_width': 12.0,
        'CBar_height': 0.50,
        'CBarLable_space': 0.25,
        'Dendrogram_width': 4.0,
        'Dendrogram_height': 12.0,
        'Dendrogram_Heatmap_gap': 0.1,
        'ScaleBar_Dendrogram_gap': 0.05,
        "linewidth_major": 1.0,
        "linewidth_minor": 0.5
    }
    ClassLabelList_x = ["Class1", "Class2", "Class3", "Class4", "Class5", "Class6"]
    ClassLabelList_y = ["Class1", "Class2", "Class3", "Class4", "Class5", "Class6"]
    TickLocList = [0, 10, 20, 30, 40, 50]

    axis_obj = construct_hmap_lines(ax_Heatmap, len_x, LineList_major, LineList_minor, hmap_params, ClassLabelList_x, ClassLabelList_y, TickLocList)

    assert isinstance(axis_obj, plt.Axes)

if __name__ == "__main__":
    pytest.main()

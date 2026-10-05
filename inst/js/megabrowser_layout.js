(function () {
  if (window.rcMegaLayoutBound) return;
  window.rcMegaLayoutBound = true;
  const pending = new WeakMap();
  function margins(workspace, plot) {
    if (!plot.closest('.mega-plot-main')) return;
    if (!plot.rcMegaMargins) plot.rcMegaMargins = {left: plot._fullLayout.margin.l, right: plot._fullLayout.margin.r};
    const small = workspace.clientWidth < 600;
    const left = small ? 30 : plot.rcMegaMargins.left;
    const right = small ? 10 : plot.rcMegaMargins.right;
    if (plot._fullLayout.margin.l !== left || plot._fullLayout.margin.r !== right) {
      window.Plotly.relayout(plot, {'margin.l': left, 'margin.r': right});
    }
  }
  function resize(workspace) {
    cancelAnimationFrame(pending.get(workspace));
    pending.set(workspace, requestAnimationFrame(function () {
      workspace.querySelectorAll('.js-plotly-plot').forEach(function (plot) {
        if (plot._fullLayout && plot.offsetWidth && plot.offsetHeight) {
          margins(workspace, plot);
          window.Plotly.Plots.resize(plot);
        }
      });
    }));
  }
  window.addEventListener('resize', function () {
    document.querySelectorAll('.mega-workspace').forEach(resize);
  });
  jQuery(document).on('shiny:value', function (event) {
    const workspace = event.target.closest('.mega-workspace');
    if (workspace && workspace.clientWidth < 600) setTimeout(function () {resize(workspace);}, 150);
  });
  document.addEventListener('input', function (event) {
    if (!event.target.matches('.mega-height-control')) return;
    const workspace = event.target.closest('.mega-workspace');
    workspace.style.setProperty('--mega-height', event.target.value + 'px');
    resize(workspace);
  });
  document.addEventListener('change', function (event) {
    if (!event.target.matches('.mega-sidebar-control,.mega-focus-control')) return;
    const workspace = event.target.closest('.mega-workspace');
    if (event.target.matches('.mega-sidebar-control')) workspace.classList.toggle('mega-hide-sidebar', !event.target.checked);
    else workspace.classList.toggle('mega-focus-heatmap', event.target.checked);
    resize(workspace);
  });
  document.addEventListener('click', function (event) {
    const button = event.target.closest('.mega-reset-view');
    if (!button) return;
    const workspace = button.closest('.mega-workspace');
    const plot = workspace.querySelector('[id$="-myPlotlyPlot"]');
    if (plot && plot._fullLayout) plot.emit('plotly_doubleclick');
  });
}());

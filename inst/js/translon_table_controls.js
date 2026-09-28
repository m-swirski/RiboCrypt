(table, checkboxId, hiddenColumns) => {
  const checkbox = document.getElementById(checkboxId);
  if (!checkbox) return;
  const update = () => {
    const visible = !checkbox.checked;
    const changed = hiddenColumns.filter(index => table.column(index).visible() !== visible);
    if (!changed.length) return;
    table.columns(changed).visible(visible, false);
    table.columns.adjust();
  };
  const event = 'change.ribocryptSimplified';
  $(checkbox).off(event).on(event, update);
  $(table.table().node()).on('destroy.dt', () => $(checkbox).off(event));
  update();
}

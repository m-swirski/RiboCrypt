(elem, _, data) => {
  let widget = HTMLWidgets.find("#" + elem.id);
  if (widget !== undefined && widget !== null) {
    let stage = widget.getStage();
    stage.tasks.onZeroOnce((x) => {
      let colorArray = data.map((hexStr) => {
        return parseInt(hexStr.replace("#", ""), 16);
      });
      let schemeId = NGL.ColormakerRegistry.addScheme(function (params) {
        this.atomColor = function (atom) {
          return colorArray[atom.resno - 1];
        };
      });
      stage.eachComponent((o) => {
        o.addRepresentation("cartoon", { color: schemeId });
        o.autoView();
      }, "structure");
    });
  };
}

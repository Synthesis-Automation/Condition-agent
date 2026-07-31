"""Reusable Qt widgets for displaying molecule and reaction drawings."""

from __future__ import annotations

from typing import Optional

from PyQt6 import QtCore, QtGui, QtSvg, QtWidgets


class StructureImageLabel(QtWidgets.QLabel):
    """Responsive label that rasterizes vector drawings at display resolution."""

    def __init__(
        self,
        *,
        placeholder: str = "Structure graph will appear here.",
        object_name: str = "structureGraph",
        minimum_height: int = 142,
    ) -> None:
        super().__init__(placeholder)
        self._placeholder = placeholder
        self._source_svg = b""
        self._source_pixmap = QtGui.QPixmap()
        self._trim_svg_white_space = True
        self._svg_max_upscale: float | None = None
        self.setObjectName(object_name)
        self.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.setMinimumWidth(360)
        self.setMinimumHeight(minimum_height)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding,
        )
        self.setStyleSheet(
            "QLabel { background: white; border: 1px solid #cbd5e0; "
            "border-radius: 4px; color: #718096; padding: 2px; }"
        )

    def set_image_bytes(
        self,
        drawing: bytes,
        *,
        trim_white_space: bool = True,
        max_upscale: float | None = None,
    ) -> bool:
        """Load image bytes with optional cropping and bounded enlargement."""
        if max_upscale is not None and max_upscale <= 0:
            raise ValueError("max_upscale must be positive")
        if b"<svg" in drawing[:512].lower():
            renderer = QtSvg.QSvgRenderer(drawing)
            if not renderer.isValid():
                self.clear_image("Unable to decode structure graph.")
                return False
            self._source_svg = bytes(drawing)
            self._source_pixmap = QtGui.QPixmap()
            self._trim_svg_white_space = trim_white_space
            self._svg_max_upscale = max_upscale
            self.setText("")
            self._refresh_pixmap()
            return True
        pixmap = QtGui.QPixmap()
        if not pixmap.loadFromData(drawing):
            self.clear_image("Unable to decode structure graph.")
            return False
        self._source_svg = b""
        self._source_pixmap = pixmap
        self._trim_svg_white_space = True
        self._svg_max_upscale = None
        self.setText("")
        self._refresh_pixmap()
        return True

    def clear_image(
        self,
        message: Optional[str] = None,
    ) -> None:
        """Remove the current drawing and show a concise placeholder."""
        self._source_svg = b""
        self._source_pixmap = QtGui.QPixmap()
        self._trim_svg_white_space = True
        self._svg_max_upscale = None
        self.setPixmap(QtGui.QPixmap())
        self.setText(message or self._placeholder)
        self.setToolTip("")

    def _refresh_pixmap(self) -> None:
        if not self._source_svg and self._source_pixmap.isNull():
            return
        available = self.contentsRect().size() - QtCore.QSize(8, 8)
        if available.width() <= 0 or available.height() <= 0:
            return
        if self._source_svg:
            self._refresh_svg(available)
            return
        self.setPixmap(
            self._source_pixmap.scaled(
                available,
                QtCore.Qt.AspectRatioMode.KeepAspectRatio,
                QtCore.Qt.TransformationMode.SmoothTransformation,
            )
        )

    def _refresh_svg(self, available: QtCore.QSize) -> None:
        renderer = QtSvg.QSvgRenderer(self._source_svg)
        source_size = renderer.defaultSize()
        if not renderer.isValid() or source_size.isEmpty():
            return
        logical_size = source_size.scaled(
            available,
            QtCore.Qt.AspectRatioMode.KeepAspectRatio,
        )
        # Extra supersampling preserves crisp strokes after whitespace trimming.
        render_scale = max(float(self.devicePixelRatioF()), 3.0)
        pixel_size = QtCore.QSize(
            max(round(logical_size.width() * render_scale), 1),
            max(round(logical_size.height() * render_scale), 1),
        )
        image = QtGui.QImage(
            pixel_size,
            QtGui.QImage.Format.Format_ARGB32_Premultiplied,
        )
        image.fill(QtCore.Qt.GlobalColor.transparent)
        painter = QtGui.QPainter(image)
        painter.setRenderHints(
            QtGui.QPainter.RenderHint.Antialiasing
            | QtGui.QPainter.RenderHint.TextAntialiasing
            | QtGui.QPainter.RenderHint.SmoothPixmapTransform
        )
        renderer.render(
            painter,
            QtCore.QRectF(
                0.0,
                0.0,
                float(pixel_size.width()),
                float(pixel_size.height()),
            ),
        )
        painter.end()
        if self._trim_svg_white_space:
            image = self._trim_white_space(
                image,
                margin=max(round(4 * render_scale), 1),
            )
        target_pixel_size = QtCore.QSize(
            max(round(available.width() * render_scale), 1),
            max(round(available.height() * render_scale), 1),
        )
        if self._svg_max_upscale is not None:
            target_pixel_size = QtCore.QSize(
                min(
                    target_pixel_size.width(),
                    max(round(image.width() * self._svg_max_upscale), 1),
                ),
                min(
                    target_pixel_size.height(),
                    max(round(image.height() * self._svg_max_upscale), 1),
                ),
            )
        image = image.scaled(
            target_pixel_size,
            QtCore.Qt.AspectRatioMode.KeepAspectRatio,
            QtCore.Qt.TransformationMode.SmoothTransformation,
        )
        pixmap = QtGui.QPixmap.fromImage(image)
        pixmap.setDevicePixelRatio(render_scale)
        self.setPixmap(pixmap)

    @staticmethod
    def _trim_white_space(
        image: QtGui.QImage,
        *,
        margin: int,
    ) -> QtGui.QImage:
        """Crop an RDKit drawing to visible content plus a small safe margin."""
        mask = image.createMaskFromColor(
            QtGui.QColor("white").rgb(),
            QtCore.Qt.MaskMode.MaskInColor,
        )
        bounds = QtGui.QRegion(QtGui.QBitmap.fromImage(mask)).boundingRect()
        if bounds.isEmpty():
            return image
        bounds.adjust(-margin, -margin, margin, margin)
        bounds = bounds.intersected(image.rect())
        return image.copy(bounds)

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        self._refresh_pixmap()


__all__ = ["StructureImageLabel"]

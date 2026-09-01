"""Focused shell presentation/routing mixin for the Qt panel."""

# ruff: noqa: F403,F405

from __future__ import annotations

from .qt_context import *  # noqa: F401,F403


class ShellMixin(QtMixinContext):
    """Cohesive GUI-only shell behavior composed into PanelController."""

    def _build_shell(self) -> None:
        self.widget.setStyleSheet(_stylesheet())
        root = self.w.QVBoxLayout(self.widget)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        header = self.w.QFrame(self.widget)
        header.setObjectName("header")
        header_layout = self.w.QHBoxLayout(header)
        header_layout.setContentsMargins(24, 18, 24, 18)
        header_layout.setSpacing(16)
        brand = self.w.QVBoxLayout()
        brand.setSpacing(2)
        brand.addWidget(_label(self.w, "STRUCTLENS / EVIDENCE BENCH", "eyebrow"))
        brand.addWidget(_label(self.w, "Integrated sequence and structure analysis", "windowTitle"))
        header_layout.addLayout(brand)
        header_layout.addStretch(1)
        self.header_status = _label(self.w, "Ready", "statusPill")
        header_layout.addWidget(self.header_status)
        self.cancel_button = _button(self.w, "Cancel", "secondaryButton")
        self.cancel_button.setVisible(False)
        self.cancel_button.clicked.connect(self._cancel_analysis)
        header_layout.addWidget(self.cancel_button)
        self.compare_button = _button(self.w, "Compare", "primaryButton")
        header_layout.addWidget(self.compare_button)
        root.addWidget(header)

        body = self.w.QWidget(self.widget)
        body_layout = self.w.QHBoxLayout(body)
        body_layout.setContentsMargins(0, 0, 0, 0)
        body_layout.setSpacing(0)
        self.sidebar = self.w.QFrame(body)
        self.sidebar.setObjectName("sidebar")
        self.sidebar.setMinimumWidth(200)
        self.sidebar.setMaximumWidth(250)
        side_layout = self.w.QVBoxLayout(self.sidebar)
        side_layout.setContentsMargins(14, 18, 14, 18)
        side_layout.setSpacing(12)
        side_layout.addWidget(_label(self.w, "WORKFLOW", "sectionKicker"))
        self.nav = self.w.QListWidget(self.sidebar)
        self.nav.setObjectName("workflowNav")
        self.nav.setSpacing(4)
        for section in SCIENTIFIC_SECTIONS:
            item = self.w.QListWidgetItem(section)
            item.setToolTip(_page_subtitle(section, standalone=self.command is None))
            self.nav.addItem(item)
        self.nav.setCurrentRow(0)
        side_layout.addWidget(self.nav, 1)
        side_layout.addWidget(
            _label(
                self.w,
                (
                    "The correspondence table is the source of truth. PyMOL only renders reversible selections."
                    if self.command is not None
                    else "The correspondence table is the source of truth. Filters and exports remain available here."
                ),
                "sidebarNote",
            )
        )
        body_layout.addWidget(self.sidebar)

        self.pages = self.w.QStackedWidget(body)
        self.pages.setObjectName("pageStack")
        body_layout.addWidget(self.pages, 1)
        root.addWidget(body, 1)

        footer = self.w.QFrame(self.widget)
        footer.setObjectName("footer")
        footer_layout = self.w.QHBoxLayout(footer)
        footer_layout.setContentsMargins(20, 8, 20, 8)
        self.footer_status = _label(self.w, "Choose two structures to begin.", "footerStatus")
        footer_layout.addWidget(self.footer_status, 1)
        self.progress = self.w.QProgressBar(footer)
        self.progress.setObjectName("analysisProgress")
        self.progress.setRange(0, 0)
        self.progress.setFixedWidth(180)
        self.progress.setVisible(False)
        footer_layout.addWidget(self.progress)
        footer_layout.addWidget(_label(self.w, "v0.3.0 · units are explicit · descriptive evidence only", "footerMeta"))
        root.addWidget(footer)

        self.compare_button.clicked.connect(self._start_analysis)

    def _add_page(self, title: str, purpose: str) -> tuple[Any, Any]:
        page = self.w.QWidget(self.pages)
        page.setObjectName(f"page{title}")
        outer = self.w.QVBoxLayout(page)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        scroll = self.w.QScrollArea(page)
        scroll.setObjectName(f"scroll{title}")
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(getattr(self.w.QFrame, "Shape", self.w.QFrame).NoFrame)
        canvas = self.w.QWidget(scroll)
        canvas.setObjectName(f"{title}PageContent")
        inner = self.w.QVBoxLayout(canvas)
        inner.setContentsMargins(30, 28, 34, 32)
        inner.setSpacing(18)
        inner.addWidget(_label(self.w, title, "pageTitle"))
        inner.addWidget(_label(self.w, purpose, "pagePurpose"))
        content = self.w.QVBoxLayout()
        content.setSpacing(16)
        inner.addLayout(content)
        inner.addStretch(1)
        scroll.setWidget(canvas)
        outer.addWidget(scroll)
        self.pages.addWidget(page)
        return page, content

    # --------------------------------------------------------------- Project

    def _wire_navigation(self) -> None:
        self.nav.currentRowChanged.connect(self.pages.setCurrentIndex)


__all__ = ["ShellMixin"]

// DashboardConsole.swift — the original tabbed engineering dashboard,
// extracted as one selectable console (F1–F6 tabs + left controls + content).

import SwiftUI

struct DashboardConsole: View {
    let supervisor: PlantSupervisor
    @Binding var activeTab: Int

    var body: some View {
        VStack(spacing: 0) {
            TabBar(activeTab: $activeTab)
                .frame(height: Theme.tabHeight)
            HStack(spacing: 0) {
                ControlsPanel(supervisor: supervisor)
                    .frame(width: Theme.controlsWidth)
                tabContent
                    .frame(maxWidth: .infinity, maxHeight: .infinity)
            }
        }
    }

    @ViewBuilder
    private var tabContent: some View {
        switch activeTab {
        case 0: OverviewTab(supervisor: supervisor)
        case 1: PrimaryTab(supervisor: supervisor)
        case 2: SecondaryTab(supervisor: supervisor)
        case 3: ReactivityTab(supervisor: supervisor)
        case 4: AlarmsTab(supervisor: supervisor)
        case 5: ICTab(supervisor: supervisor)
        default: EmptyView()
        }
    }
}
